import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

import crn_core


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


def run(args):
    assert args.nrobots > 1, 'There must be at least 2 robot species.'

    if not CreateDirectory(args.directory, args.force):
        return
    print ''

    # Setup CRN (do not add the collaboration state when computing CME).
    #crn = crn_core.TaskCRNSetup(args.nrobots, args.stochkit_path)
    crn = crn_core.TaskCRNSetup(args.nrobots, args.stochkit_path, alpha=args.alpha, beta=args.beta)
    
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

    # Show plots if requested.
    if args.show_plots:
        plotter = crn_core.Plotter(output)
        # Plot average.
        #plotter.AverageTrajectory()
        # Plot distribution of the species populations.
        plotter.Distributions_Single(from_timestamp=args.duration / 2.0)

    # Show plots of the observable data.
    if args.show_observable_plots:
        aggregated_output = crn_core.TaskCRNExtractObservableData(output)
        plotter = crn_core.Plotter(aggregated_output)
        # Plot average.
        #plotter.AverageTrajectory()
        # Plot distribution of the species populations.
        plotter.Distributions(from_timestamp=args.duration / 2.0)

    if args.show_plots or args.show_observable_plots:
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
    parser.add_argument('--cme', action='store_true', help='If set, uses CMEPy instead of StochKit')
    
    parser.add_argument('--alpha', metavar='RATE', type=float, nargs='+', required=True, help='Alpha parameters (for all species)')
    parser.add_argument('--beta', metavar='RATE', type=float, nargs='+', required=True, help='Beta parameters (for all species)')
    
    # Run.
    run(parser.parse_args())
