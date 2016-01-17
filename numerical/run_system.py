import argparse
import cmepy as cme
import cmepy.recorder
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

from crn import CRN
from plot_utils import Plotter

# Default rate for all species for EXPLORE -> WAIT.
_DEFAULT_ALPHA = 0.2

# Default rate for all species for WAIT -> EXPLORE (without collaboration).
_DEFAULT_BETA = 0.1


def CreateDirectory(directory, force=False):
    if force and os.path.exists(directory):
        print 'Deleting previous directory %s.' % directory
        shutil.rmtree(directory)
    print 'Preparing directory %s.' % directory
    try:
        os.makedirs(directory)
    except os.error:
        print 'Cannot create directory %s. Make sure it does not exist.' % directory
        return False
    return True


def SetupCRN(nrobots, stochkit_binary, add_collaboration_state=True,
             alpha=_DEFAULT_ALPHA, beta=_DEFAULT_BETA):
    crn = CRN(stochkit_binary=stochkit_binary)
    nspecies = len(nrobots)
    N = [str(i) for i in xrange(1, nspecies + 1)]  # Species names (i.e., 1, 2, 3, ...).
    # Species go from EXPLORE to WAIT.
    # We will use 'e' for explore and 'w' for wait.
    # A condition like '124w' means that one robot of species 1, 2 and 4 are waiting.
    # TODO: Use a configured alpha and beta for each species.
    waiting_conditions = set()
    for i in xrange(nspecies):
        crn.AddReaction([N[i] + 'e'], [N[i] + 'w'], alpha)
        crn.AddReaction([N[i] + 'w'], [N[i] + 'e'], beta)
        waiting_conditions.add((i,))
    # Robots in the waiting state can be joined by new robots.
    while waiting_conditions:
        new_waiting_conditions = set()
        for condition in waiting_conditions:
            condition_string = ''.join([N[j] for j in condition]) + 'w'
            for i in xrange(nspecies):
                # A robot cannot join a condition where its species is already present.
                if i in condition:
                    continue
                # Add robot to the condition (and sort robot types).
                new_condition = tuple(sorted(condition + (i,)))
                # Special case when all species are present.
                if len(new_condition) == nspecies:
                    product = [(N[j] + 'e') for j in xrange(nspecies)]
                    if add_collaboration_state:
                        product += ['c']
                    crn.AddReaction([condition_string, N[i] + 'e'], product, alpha)
                else:
                    new_condition_string = ''.join([N[j] for j in new_condition]) + 'w'
                    # Add forward and backward reactions
                    crn.AddReaction([condition_string, N[i] + 'e'], [new_condition_string], alpha)
                    crn.AddReaction([new_condition_string], [condition_string, N[i] + 'e'], beta)
                    new_waiting_conditions.add(new_condition)
        waiting_conditions = new_waiting_conditions
    # Set population of robots (starting in exploration).
    for i, n in enumerate(nrobots):
        crn.SetInitialPopulation(N[i] + 'e', n)
    return crn


# This function is very hacky. In particular it assumes that the conditions are
# named like: 'Xe', 'Xw', 'XYw', ... and thus uses the length and the last letter
# of the condition to aggregate the data.
def ExtractObservableData(output_from_run):
    # Helper function for sum relevant conditions in each aggregate.
    def _ConditionsCountsFunction(aggregate):
        return lambda *x: reduce(lambda v, w: v + x[w], aggregate, 0)

    if len(output_from_run) == 3:
        # StochKit data.
        conditions, timestamps, data = output_from_run
        n = max([len(c) for c in conditions]) - 1  # Hacky: Get largest aggregate size from the condition names.
        # Aggregated data for 'Xe', 'Xw', 'XYw', ...
        aggregated_data = np.zeros((data.shape[0], data.shape[1], n + 1))
        for i, c in enumerate(conditions):
            if c == 'c':
                continue
            if c.endswith('e'):
                aggregated_data[:, :, 0] += data[:, :, i]
                continue
            aggregated_data[:, :, len(c) - 1] += data[:, :, i]
        aggregated_conditions = ['?e'] + [('?' * i + 'w') for i in xrange(1, n + 1)]
        return aggregated_conditions, timestamps, aggregated_data
    elif len(output_from_run) == 2:
        # CME data.
        conditions, recorder = output_from_run
        # List all observable conditions.
        n = max([len(c) for c in conditions]) - 1  # Hacky: Get largest aggregate size from the condition names.
        aggregates = [[] for _ in xrange(n + 1)]
        for i, c in enumerate(conditions):
            if c == 'c':
                continue
            if c.endswith('e'):
                aggregates[0].append(i)
                continue
            aggregates[len(c) - 1].append(i)
        aggregated_conditions = ['?e'] + [('?' * i + 'w') for i in xrange(1, n + 1)]
        new_recorder = cmepy.recorder.create((aggregated_conditions, [_ConditionsCountsFunction(a) for a in aggregates]))
        conditions_tuple = tuple(conditions)
        measurement = recorder[conditions_tuple]
        time_steps = measurement.times
        for i, t in enumerate(time_steps):
            new_recorder.write(t, measurement.distributions[i])
        return aggregated_conditions, new_recorder
    return None


# This function extracts the performance data.
def ExtractPerformanceData(output_from_run):
    assert len(output_from_run) == 3, 'ExtractPerformanceData() can only be used on StochKit output.'
    conditions, timestamps, data = output_from_run
    aggregated_data = np.zeros((data.shape[0], data.shape[1], 1))
    for i, c in enumerate(conditions):
        if c == 'c':
            aggregated_data[:, :, 0] += data[:, :, i]
    aggregated_conditions = ['c']
    return aggregated_conditions, timestamps, aggregated_data


def run(args):
    assert args.nrobots > 1, 'There must be at least 2 robot species.'

    if not CreateDirectory(args.directory, args.force):
        return
    print ''

    # Setup CRN (do not add the collaboration state when computing CME).
    crn = SetupCRN(args.nrobots, args.stochkit_path, add_collaboration_state=not args.cme)
    print 'Running simulation on:'
    print crn, '\n'

    # Run simulation.
    if args.cme:
        # Try to be more precise at the beginning.
        time_steps = np.concatenate((
            np.arange(0., args.duration / 4.0, 0.1),             # Every 0.1 seconds
            np.arange(args.duration / 4.0, args.duration, 0.5),  # Every 0.5 seconds.
            np.array([args.duration])  # Last timestep.
        ))
        output = crn.CME(time_steps)
    else:
        output = crn.Simulate(duration=args.duration, nruns=args.nruns, output_directory=args.directory)
    # One can compute the observable data distribution as follows:
    # aggregated_output = ExtractObservableData(output)
    # observable_distribution = crn.BuildDistribution(from_timestamp=args.duration / 2.0)
    #
    # The observable distribution can then be compared over two runs of the crn:
    # crn.CompareDistributions(observable_distribution_a, observable_distribution_b)
    #
    # The performance can also be extracted as:
    # performance_output = ExtractPerformanceData(output)

    # Show plots if requested.
    if args.show_plots:
        plotter = Plotter(output, exclude_list=('c',))  # Do not show collaborations in the plots.
        # Plot average.
        plotter.AverageTrajectory()
        # Plot distribution of the species populations.
        plotter.Distributions(from_timestamp=args.duration / 2.0)

    # Show plots of the observable data.
    if args.show_observable_plots:
        aggregated_output = ExtractObservableData(output)
        plotter = Plotter(aggregated_output)
        # Plot average.
        plotter.AverageTrajectory()
        # Plot distribution of the species populations.
        plotter.Distributions(from_timestamp=args.duration / 2.0)

    # Show plots of the performance.
    if args.show_performance_plots and not args.cme:
        aggregated_output = ExtractPerformanceData(output)
        plotter = Plotter(aggregated_output)
        # Plot average.
        plotter.AverageTrajectory()

    if args.show_plots or args.show_observable_plots or args.show_performance_plots:
        plt.show()


if __name__ == '__main__':
    # Setup arguments.
    parser = argparse.ArgumentParser(description='Runs StochKit on the task allocation problem')
    parser.add_argument('--nrobots', metavar='N', type=int, nargs='+', required=True, help='Number of robots in each species')
    parser.add_argument('--directory', metavar='PATH', action='store', default='/tmp/run_system_output', help='Path where the data is stored')
    parser.add_argument('--stochkit_path', metavar='PATH', action='store', default='./StochKit2.0.11/ssa', help='Path where the stochkit binary is stored')
    parser.add_argument('--nruns', metavar='N', type=int, action='store', default=100, help='Number of simulation runs')
    parser.add_argument('--duration', metavar='SECONDS', type=float, action='store', default=20., help='Duration in seconds of each simulation run')
    parser.add_argument('--force', action='store_true', help='If set, the directory is overwritten')
    parser.add_argument('--show_plots', action='store_true', help='If set, plots are shown after the simulations')
    parser.add_argument('--show_observable_plots', action='store_true', help='If set, plots are shown after the simulations with observable data aggregated')
    parser.add_argument('--show_performance_plots', action='store_true', help='If set, plots are shown after the simulations with performance data aggregated')
    parser.add_argument('--cme', action='store_true', help='If set, uses CMEPy instead of StochKit')
    # Run.
    run(parser.parse_args())
