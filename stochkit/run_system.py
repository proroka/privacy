import argparse
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


def SetupCRN(nrobots, stochkit_binary):
    crn = CRN(stochkit_binary=stochkit_binary)
    nspecies = len(nrobots)
    N = [str(i) for i in xrange(1, nspecies + 1)]  # Species names (i.e., 1, 2, 3, ...).
    # Species go from EXPLORE to WAIT.
    # We will use 'e' for explore and 'w' for wait.
    # A condition like '124w' means that one robot of species 1, 2 and 4 are waiting.
    # TODO: Use a configured alpha and beta for each species.
    waiting_conditions = set()
    for i in xrange(nspecies):
        crn.AddReaction([N[i] + 'e'], [N[i] + 'w'], _DEFAULT_ALPHA)
        crn.AddReaction([N[i] + 'w'], [N[i] + 'e'], _DEFAULT_BETA)
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
                    crn.AddReaction([condition_string, N[i] + 'e'], [(N[j] + 'e') for j in xrange(nspecies)] + ['c'], _DEFAULT_ALPHA)
                else:
                    new_condition_string = ''.join([N[j] for j in new_condition]) + 'w'
                    # Add forward and backward reactions
                    crn.AddReaction([condition_string, N[i] + 'e'], [new_condition_string], _DEFAULT_ALPHA)
                    crn.AddReaction([new_condition_string], [condition_string, N[i] + 'e'], _DEFAULT_BETA)
                    new_waiting_conditions.add(new_condition)
        waiting_conditions = new_waiting_conditions
    # Set population of robots (starting in exploration).
    for i, n in enumerate(nrobots):
        crn.SetInitialPopulation(N[i] + 'e', n)
    return crn


def run(args):
    assert args.nrobots > 1, 'There must be at least 2 robot species.'

    if not CreateDirectory(args.directory, args.force):
        return
    print ''

    # Setup CRN.
    crn = SetupCRN(args.nrobots, args.stochkit_path)
    print 'Running simulation on:'
    print crn, '\n'

    # Run simulation.
    output = crn.Run(duration=args.duration, nruns=args.nruns, output_directory=args.directory)

    # Show plots if requested.
    if args.show_plots:
        plotter = Plotter(output, exclude_list=('c',))  # Do not show collaborations in the plots.
        # Plot average.
        plotter.AverageTrajectory()
        # Plot distribution of the species populations.
        plotter.Distributions(from_timestamp=args.duration / 2.0)

    # Show plots of the observable data.
    if args.show_observable_plots:
        conditions, timestamps, data = output
        n = max([len(c) for c in conditions]) - 1  # Hacky: Get largest aggregate size from the condition names.
        # Aggregated data for 'Xe', 'Xw', 'XYw', ...
        aggregated_data = np.zeros((args.nruns, timestamps.shape[0], n + 1))
        for i, c in enumerate(conditions):
            if c == 'c':
                continue
            if c.endswith('e'):
                aggregated_data[:, :, 0] += data[:, :, i]
                continue
            aggregated_data[:, :, len(c) - 1] += data[:, :, i]
        aggregated_conditions = ['Exploring'] + [('Tasks with %d robots' % i) for i in xrange(1, n + 1)]
        plotter = Plotter((aggregated_conditions, timestamps, aggregated_data))
        # Plot average.
        plotter.AverageTrajectory()
        # Plot distribution of the species populations.
        plotter.Distributions(from_timestamp=args.duration / 2.0)

    if args.show_plots or args.show_observable_plots:
        plt.show()


if __name__ == '__main__':
    # Setup arguments.
    parser = argparse.ArgumentParser(description='Runs StochKit on the task allocation problem.')
    parser.add_argument('--nrobots', metavar='N', type=int, nargs='+', required=True, help='Number of robots in each species')
    parser.add_argument('--directory', metavar='PATH', action='store', default='/tmp/run_system_output', help='Path where the data is stored')
    parser.add_argument('--stochkit_path', metavar='PATH', action='store', default='./StochKit2.0.11/ssa', help='Path where the stochkit binary is stored')
    parser.add_argument('--nruns', metavar='N', type=int, action='store', default=100, help='Number of simulation runs.')
    parser.add_argument('--duration', metavar='SECONDS', type=float, action='store', default=20., help='Duration in seconds of each simulation run.')
    parser.add_argument('--force', action='store_true', help='If set, the directory is overwritten')
    parser.add_argument('--show_plots', action='store_true', help='If set, plots are shown after the simulations')
    parser.add_argument('--show_observable_plots', action='store_true', help='If set, plots are shown after the simulations with observable data aggregated')
    # Run.
    run(parser.parse_args())
