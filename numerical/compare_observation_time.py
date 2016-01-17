import argparse
import matplotlib.pyplot as plt
import numpy as np
import pickle

import crn_core

# Default rate for all species for EXPLORE -> WAIT.
_DEFAULT_ALPHA = 0.2

# Default rate for all species for WAIT -> EXPLORE (without collaboration).
_DEFAULT_BETA = 0.1

# Default duration
_DEFAULT_DURATION = 5.

# Number of points per axis in the graph.
_NPOINTS = 10

# To only show the first call to RunSystem.
first_run = True


def RunSystem(populations, taus, args):
    global first_run

    # Ignore collaboration rate.
    crn_model = crn_core.TaskCRNSetup(populations, None,  alpha=args.alpha, beta=args.beta)
    # Run simulation.
    time_steps = np.array(taus)
    output = crn_model.CME(time_steps)
    aggregated_output = crn_core.TaskCRNExtractObservableData(output)
    if args.plot_first_run and first_run:
        print 'Showing plot for the population:', populations
        p = crn_core.Plotter(aggregated_output)
        p.AverageTrajectory()
        p.Distributions(from_timestamp=taus[-1] * 0.95)
        plt.show()
        first_run = False
    return crn_core.BuildDistributions(aggregated_output)


def GetAllEpsilons(args):
    # Base population (all species have the same number of robots).
    base_populations = np.array(args.nrobots)
    nspecies = len(args.nrobots)
    # Alternative databases. For each species, we can remove a robots or switch its species.
    alternative_population_offsets = []
    for s in xrange(nspecies):
        # Remove robot.
        remove_offset_array = np.zeros((nspecies,)).astype(int)
        remove_offset_array[s] = -1
        # Removing a robot makes no sense as the number of robots is observed.
        # i.e., alternative_population_offsets.append(remove_offset_array)
        # Swich to another team.
        for other_s in xrange(nspecies):
            if other_s == s:
                continue
            add_offset_array = np.zeros((nspecies,)).astype(int)
            add_offset_array[other_s] = 1.
            alternative_population_offsets.append(remove_offset_array + add_offset_array)
    taus = np.linspace(0., args.duration, args.npoints).tolist()
    # Keeps track of all taus.
    min_epsilons = {}
    print 'Analyzing taus =', taus
    base_distribution = RunSystem(base_populations, taus, args)
    # Alternative database.
    max_min_epsilons = [0.] * len(taus)
    for alternative_populations in alternative_population_offsets:
        alternative_distribution = RunSystem(alternative_populations + base_populations, taus, args)
        for i, tau in enumerate(taus):
            # Compute difference.
            min_epsilon = crn_core.CompareDistributions(base_distribution[i], alternative_distribution[i], prune=1e-3)
            # Store difference.
            max_min_epsilons[i] = max(max_min_epsilons[i], min_epsilon)
    for i, tau in enumerate(taus):
        min_epsilons[tau] = max_min_epsilons[i]
    return min_epsilons


def run(args):
    if args.load_epsilons:
        print 'Loading precomputed epsilons...'
        min_epsilons = pickle.load(open(args.load_epsilons, 'r'))
    else:
        min_epsilons = GetAllEpsilons(args)
        if args.save_epsilons:
            print 'Saving epsilons...'
            pickle.dump(min_epsilons, open(args.save_epsilons, 'w'))
    # Plot it!
    values = np.array(sorted(min_epsilons.iteritems()))
    plt.figure()
    plt.plot(values[:, 0], values[:, 1], linewidth=2, color='b', label='\epsilon')
    plt.title('N = %s, alpha = %s, beta = %s' % (str(args.nrobots), str(args.alpha), str(args.beta)))
    plt.xlabel('tau')
    plt.ylabel('Minimum acheivable \epsilon')
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Outputs minimum epsilon for differential privacy')
    parser.add_argument('--nrobots', metavar='N', type=int, nargs='+', required=True, help='Number of robots in each species')
    parser.add_argument('--alpha', metavar='RATE', type=float, nargs='+', required=True, help='Alpha parameters (for all species)')
    parser.add_argument('--beta', metavar='RATE', type=float, nargs='+', required=True,  help='Beta parameters (for all species)')
    parser.add_argument('--duration', metavar='SECONDS', type=float, action='store', default=_DEFAULT_DURATION, help='Run duration')
    parser.add_argument('--nruns', metavar='N', type=int, action='store', default=1000, help='Number of StochKit runs')
    parser.add_argument('--npoints', metavar='N', type=int, action='store', default=_NPOINTS, help='Number datapoints')
    parser.add_argument('--save_epsilons', metavar='FILE', type=str, action='store', help='If set, saves the min-epsilons in the specified file')
    parser.add_argument('--load_epsilons', metavar='FILE', type=str, action='store', help='If set, load the min-epsilons from the specified file')
    parser.add_argument('--plot_first_run', action='store_true', help='If set, plots the first run')

    # Run.
    run(parser.parse_args())
