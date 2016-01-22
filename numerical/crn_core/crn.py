import argparse
import cmepy as cme
import cmepy.model
import cmepy.recorder
import cmepy.fsp.solver
import cmepy.fsp.support_expander
import cmepy.domain
import cmepy.statistics
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import subprocess
import tempfile
from xml.dom import minidom
from xml.etree import ElementTree
from xml.etree.ElementTree import Element, SubElement, Comment

import plot_utils


# Name of the StochKit XML file.
_STOCHKIT_FILENAME = 'stochkit.xml'
# Name of folders created by StochKit.
_STOCHKIT_STATS_OUTPUT_MEAN = 'stochkit_output/stats/means.txt'
_STOCHKIT_STATS_OUTPUT_STD = 'stochkit_output/stats/variances.txt'
_STOCHKIT_TRAJECTORIES_OUTPUT = 'stochkit_output/trajectories/trajectory%d.txt'

# Probability value below which distribution difference is not computed due
# to numerical imprecision.
_NU = 1e-5

# Defines a complex.
class Complex(object):
    def __init__(self, species):
        # Build a dictionary of <species, counts>.
        self.counts = {}
        for s in species:
            if s not in self.counts:
                # We define the CRN species as normal strings.
                self.counts[s] = 0
            self.counts[s] += 1

    def __str__(self):
        if self.counts:
            return ' + '.join((('' if c == 1 else str(c)) + s) for s, c in self.counts.iteritems())
        return '(null)'

    def AddXML(self, node, node_type):
        # Type should be either 'Reactants' or 'Products'.
        assert node_type == 'Reactants' or node_type == 'Products'
        subnode = SubElement(node, node_type)
        for s, c in self.counts.iteritems():
            SubElement(subnode, 'SpeciesReference', {'id': s, 'stoichiometry': '%d' % c})
        if not self.counts:
            subnode.append(Comment('no ' + node_type.lower()))


# Defines a reaction.
class Reaction(object):
    def __init__(self, inputs, outputs, rate):
        self.inputs = inputs
        self.outputs = outputs
        self.rate = rate

    def __str__(self):
        return str(self.inputs) + ' -> ' + str(self.outputs)

    def AddXML(self, node, id, ):
        reaction_node = SubElement(node, 'Reaction')
        SubElement(reaction_node, 'Id').text = 'R' + str(id + 1)
        SubElement(reaction_node, 'Description').text = self.__str__()
        SubElement(reaction_node, 'Type').text = 'mass-action'
        SubElement(reaction_node, 'Rate').text = '%g' % self.rate
        self.inputs.AddXML(reaction_node, 'Reactants')
        self.outputs.AddXML(reaction_node, 'Products')


# Defines a CRN.
class CRN(object):
    def __init__(self, stochkit_binary='ssa'):
        self.species = set()  # Set of strings.
        self.species_descriptions = {}  # Species can have a description (defaults to 'Species X').
        self.species_population = {}  # Species have an initial population (defaults to 0).
        self.reactions = []   # List of Reactions.
        self.stochkit_binary = stochkit_binary
        self.ordered_species = []  # Always keep a consistent species order.

    def AddReaction(self, inputs, outputs, rate):
        self.species |= set(inputs)
        self.species |= set(outputs)
        self.ordered_species = sorted(self.species)
        self.reactions.append(Reaction(Complex(inputs), Complex(outputs), rate))

    def SetSpeciesDescription(self, species, description):
        assert species in self.species, 'No reactions with species "%s". Description cannot be set.' % species
        self.species_descriptions[species] = description

    def SetInitialPopulation(self, species, population):
        assert species in self.species, 'No reactions with species "%s". Population cannot be set.' % species
        self.species_population[species] = population

    def __str__(self):
        return '\n'.join(str(r) for r in self.reactions)

    def _XMLString(self, description):
        # This function creates a StochKit compatible XML code.
        model_node = Element('Model')
        model_node.append(Comment('Generated by crn.py'))
        SubElement(model_node, 'Description').text = description
        # General info.
        SubElement(model_node, 'NumberOfReactions').text = str(len(self.reactions))
        SubElement(model_node, 'NumberOfSpecies').text = str(len(self.species))
        # Reactions.
        reactionslist_node = SubElement(model_node, 'ReactionsList')
        for i, reaction in enumerate(self.reactions):
            reaction.AddXML(reactionslist_node, i)
        # Species.
        specieslist_node = SubElement(model_node, 'SpeciesList')
        for i, species in enumerate(self.ordered_species):
            species_node = SubElement(specieslist_node, 'Species')
            SubElement(species_node, 'Id').text = species
            SubElement(species_node, 'Description').text = (
                self.species_descriptions[species] if species in self.species_descriptions else
                ('Species %s' % species))
            SubElement(species_node, 'InitialPopulation').text = (
                '%d' % (self.species_population[species] if species in self.species_population else 0))
        # Small trick to make the XML pretty.
        raw_output = ElementTree.tostring(model_node, 'utf-8')
        return minidom.parseString(raw_output).toprettyxml(indent="  ")

    # Runs the complete simulation using StochKit and return 3 values:
    # (a) The list of species
    # (b) The timestamps at which the simulation is snapshot
    # (c) The species population at these timestamps as a 3D tensor with dimensions nruns x ndatapoints x nspecies.
    def Simulate(self, duration=10.0, nruns=1, ndatapoints=None, output_directory=None):
        # Create XML file (overwrite old file if any).
        # If output_directory is None, create a temporary directory.
        cleandir = False
        if output_directory is None:
            output_directory = tempfile.mkdtemp()
            cleandir = True
        xmlpath = os.path.join(output_directory, _STOCHKIT_FILENAME)
        with open(xmlpath, 'w') as fp:
            fp.write(self._XMLString('Autogenerated CRN'))
        # Keep one datapoint per time unit if ndatapoints is not set.
        if ndatapoints is None:
            ndatapoints = int(duration) + 1
        # Prepare command.
        commandline = [
            self.stochkit_binary,
            '-m', xmlpath,
            '-t', '%g' % duration,
            '-r', '%d' % nruns,
            '-i', '%d' % (ndatapoints - 1),
            '--keep-trajectories',
            '-f']
        print 'Executing:', ' '.join(commandline)
        if subprocess.call(commandline) != 0:
            raise RuntimeError('Error occur while running StochKit')
        # Gather trajectories back into numpy arrays.
        timestamps = np.empty((ndatapoints, 1))
        data = np.empty((nruns, ndatapoints, len(self.species)))
        for run in xrange(nruns):
            run_filename = os.path.join(output_directory, _STOCHKIT_TRAJECTORIES_OUTPUT % run)
            d = np.loadtxt(run_filename)
            data[run, :, :] = d[:, 1:]
            if run == 0:
                timestamps = d[:, 0]
        if cleandir:
            shutil.rmtree(output_directory)
        return self.ordered_species, timestamps, data

    def _CMEModel(self, description):
        # Utilities lambda generators.
        def _SpeciesCountsFunction(i):
            return lambda *x: x[i]

        def _PropensitiesFunction(reaction, species_counts, species_to_index):
            def action_mass(*x):
                r = reaction.rate
                for s, c in reaction.inputs.counts.iteritems():
                    r *= species_counts[species_to_index[s]](*x) ** float(c)
                return r
            return action_mass

        def _GetTransitions(reaction, species_to_index):
            transitions = [0]*len(species_to_index)
            for s, c in reaction.inputs.counts.iteritems():
                transitions[species_to_index[s]] -= c
            for s, c in reaction.outputs.counts.iteritems():
                transitions[species_to_index[s]] += c
            return transitions

        species_to_index = dict((s, i) for i, s in enumerate(self.ordered_species))
        species_counts = [_SpeciesCountsFunction(i) for i in xrange(len(self.ordered_species))]
        reactions = []
        propensities = []
        transitions = []
        for i, reaction in enumerate(self.reactions):
            reactions.append(str(reaction))
            propensities.append(_PropensitiesFunction(reaction, species_counts, species_to_index))
            transitions.append(_GetTransitions(reaction, species_to_index))
        initial_state = [0]*len(species_to_index)
        for s, c in self.species_population.iteritems():
            initial_state[species_to_index[s]] = c
        initial_state = tuple(initial_state)
        return cme.model.create(name=description, species=self.ordered_species, species_counts=species_counts,
                                reactions=reactions, propensities=propensities, transitions=transitions,
                                initial_state=initial_state)

    # Runs the complete simulation using CMEPy and return 4 values:
    # (a) The list of species
    # (b) The cmepy.recorder object. You can use BuildDistribution and the Plotter utility on it.
    # See http://fcostin.github.io/cmepy/modules/recorder.html#module-recorder for more details on it.
    def CME(self, time_steps, distribution_precision=1e-2, verbose=False):
        model = self._CMEModel('Autogenerated CRN')
        initial_states = cme.domain.from_iter((model.initial_state, ))
        # Create expander for FSP expansion strategy.
        expander = cme.fsp.support_expander.SupportExpander(
            model.transitions, depth=1, epsilon=1.e-5)
        # Create FSP solver/
        fsp_solver = cme.fsp.solver.create(model, initial_states, expander)
        # Error of the solution at the final time to be bounded.
        num_steps = np.size(time_steps)
        max_error_per_step = distribution_precision / float(num_steps)
        # Create recorder to record species counts.
        recorder = cmepy.recorder.create((model.species, model.species_counts))
        print 'Executing: CME'
        for i, t in enumerate(time_steps):
            if i % max(len(time_steps) / 20, 1) == 0 and verbose:
                print '%d%%' % (i * 100 / len(time_steps))
            fsp_solver.step(t, max_error_per_step)
            # Record the solution at this timestep.
            recorder.write(t, fsp_solver.y[0])
        print 'done!'
        return self.ordered_species, recorder


# This function builds an N-D tensor (stored as a sparse numpy matrix)
# That matches the distribution and axis order. This matrix can be
# used to compare two steady-state distributions.
def BuildDistribution(output_from_run, from_timestamp=0.):
    if len(output_from_run) == 3:
        # StochKit.
        counts = {}
        species, timestamps, trajectories = output_from_run
        trajectories = trajectories[:, timestamps >= from_timestamp, :]
        nruns = trajectories.shape[0]
        ntimestamps = trajectories.shape[1]
        total_counts = 0
        for i in xrange(nruns):
            for j in xrange(ntimestamps):
                indices = tuple(trajectories[i, j, :].astype(int).tolist())
                if indices not in counts:
                    counts[indices] = 0
                counts[indices] += 1
                total_counts += 1
        for k, v in counts.iteritems():
            counts[k] = float(counts[k]) / float(total_counts)
        return counts
    elif len(output_from_run) == 2:
        # CMEPy.
        species, recorder = output_from_run
        return recorder[tuple(species)].distributions[-1]
    raise ValueError('Unable to detect the type of data.')


# Same as above but outputs the distributions of all timesteps.
def BuildDistributions(output_from_run):
    if len(output_from_run) == 3:
        raise ValueError('Unable to process StochKit data. Use CMEPy.')
    elif len(output_from_run) == 2:
        # CMEPy.
        species, recorder = output_from_run
        return recorder[tuple(species)].distributions
    raise ValueError('Unable to detect the type of data.')


# Compares two distributions computed by BuildDistribution().
def CompareDistributions(A, B, smooth=_NU):
1e-5ll_values = []
   for k, v in A.iteritems():
       if k in B:
           all_values.append(np.abs(np.log(A[k] + smooth) - np.log(B[k] + smooth)))
       elif k not in B:
           all_values.append(np.abs(np.log(A[k] + smooth) - np.log(smooth)))
       # Ignore when one is larger and the other one is smaller or when both are smaller.
   for k, v in B.iteritems():
       if k not in A:
           all_values.append(np.abs(np.log(smooth) - np.log(B[k] + smooth)))
   return max(all_values)


if __name__ == '__main__':
    # Parsing arguments to get the path to the StochKit binary.
    parser = argparse.ArgumentParser(description='Runs StochKit on a test CRN.')
    parser.add_argument('--stochkit_path', metavar='PATH', action='store', default='./StochKit2.0.11/ssa', help='Path where the stochkit binary is stored')
    args = parser.parse_args()
    # Testing.
    crn = CRN(stochkit_binary=args.stochkit_path)
    crn.AddReaction(['A', 'B'], ['C'], 1.)  # A + B -> C
    crn.AddReaction(['C'] * 2, ['D'], 0.5)  # 2C -> D
    crn.AddReaction(['D'], [], 2.)          # D -> (nothing)
    crn.AddReaction([], ['A', 'B'], 1.)     # (nothing) -> A + B
    crn.SetInitialPopulation('A', 10)
    print 'Running simulation on:'
    print crn, '\n'
    # Run simulation.
    output = crn.Simulate(duration=10., nruns=1000)
    # Plot average.
    plotter = plot_utils.Plotter(output)
    plotter.AverageTrajectory()
    # Plot distribution of the species populations after 5 seconds.
    plotter.Distributions(from_timestamp=5.)
    # Re-run StochKit to test the distribution comparison function.
    d1 = BuildDistribution(output, from_timestamp=5.)
    output = crn.Simulate(duration=50., nruns=100)
    d2 = BuildDistribution(output, from_timestamp=5.)
    print 'Comparison of two runs of StochKit:', CompareDistributions(d1, d2)

    # CME (we give more precision at the beginning of the simulation).
    time_steps = np.concatenate((
        np.linspace(0., 2., 21),  # Every 0.1 seconds
        np.linspace(3., 10., 8)  # Every seconds.
    ))
    output = crn.CME(time_steps)
    # Plot average.
    plotter = plot_utils.Plotter(output)
    plotter.AverageTrajectory()
    # Plot distribution of the species populations at the end.
    plotter.Distributions()
    # Compare with StochKit.
    d3 = BuildDistribution(output, from_timestamp=5.)
    print 'Comparison of CME with StochKit:', CompareDistributions(d1, d3)
    output = crn.CME(time_steps)
    d4 = BuildDistribution(output, from_timestamp=5.)
    print 'Comparison of two CME runs:', CompareDistributions(d3, d4, prune=1e-6)
    # Show the plots.
    plt.show()
