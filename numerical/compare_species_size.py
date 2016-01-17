import argparse
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pickle
from scipy import interpolate

import crn
import plot_utils
import run_system as system

# Default rate for all species for EXPLORE -> WAIT.
_DEFAULT_ALPHA = 0.2

# Default rate for all species for WAIT -> EXPLORE (without collaboration).
_DEFAULT_BETA = 0.1

# Default duration
_DEFAULT_DURATION = 5.

# To only show the first call to RunSystem.
first_run = True


def RunSystem(populations, args):
    global first_run

    # Ignore collaboration rate.
    crn_model = system.SetupCRN(populations, args.stochkit_path, add_collaboration_state=False,
                                alpha=args.alpha, beta=args.beta)
    # Run simulation.
    if args.cme:
        # Try to be more precise at the beginning.
        time_steps = np.concatenate((
            np.arange(0., args.duration / 4., .1),             # Every 0.1 seconds
            np.arange(args.duration / 4., args.duration, .5),  # Every 0.5 seconds.
            np.array([args.duration])  # Last timestep.
        ))
        output = crn_model.CME(time_steps)
    else:
        output = crn_model.Simulate(duration=args.duration, nruns=args.nruns)
    aggregated_output = system.ExtractObservableData(output)
    if args.plot_first_run and first_run:
        print 'Showing plot for the population:', populations
        p = plot_utils.Plotter(aggregated_output)
        p.AverageTrajectory()
        p.Distributions(from_timestamp=args.duration / 2.0)
        plt.show()
        first_run = False
    return crn.BuildDistribution(aggregated_output, from_timestamp=args.duration / 2.0)


def GetAllEpsilons(args):
    # 1 species is fixed at args.nrobots, while the others vary in numbers.
    populations_to_test = itertools.product(xrange(1, args.nrobots + 1), repeat=args.nspecies - 1)
    populations_to_test = [np.array((args.nrobots,) + p).astype(int) for p in populations_to_test]
    # Alternative databases. For each species, we can remove a robots or switch its species.
    alternative_population_offsets = []
    for s in xrange(args.nspecies):
        # Remove robot.
        remove_offset_array = np.zeros((args.nspecies,)).astype(int)
        remove_offset_array[s] = -1
        # Removing a robot makes no sense as the number of robots is observed.
        # i.e., alternative_population_offsets.append(remove_offset_array)
        # Swich to another team.
        for other_s in xrange(args.nspecies):
            if other_s == s:
                continue
            add_offset_array = np.zeros((args.nspecies,)).astype(int)
            add_offset_array[other_s] = 1.
            alternative_population_offsets.append(remove_offset_array + add_offset_array)
    # Keeps track of all distributions.
    # Since the system is symmetric, we store in sorted population numbers.
    # So that we do not test both (N_1, N_2, N_3) and (N_1, N_3, N_2).
    distributions_cache = {}
    min_epsilons_cache = {}
    min_epsilons = {}
    for populations in populations_to_test:
        unsorted_base_populations_tuple = tuple(populations.tolist())
        base_populations_tuple = tuple(sorted(populations.tolist()))
        # We might have computed this population before.
        if base_populations_tuple in distributions_cache:
            base_distribution = distributions_cache[base_populations_tuple]
        else:
            base_distribution = RunSystem(base_populations_tuple, args)
            distributions_cache[base_populations_tuple] = base_distribution
        # Alternative database.
        max_min_epsilon = 0.
        for alternative_populations in alternative_population_offsets:
            alternative_populations_tuple = tuple(sorted((alternative_populations + populations).tolist()))
            if alternative_populations_tuple in distributions_cache:
                alternative_distribution = distributions_cache[alternative_populations_tuple]
            else:
                alternative_distribution = RunSystem(alternative_populations_tuple, args)
                distributions_cache[alternative_populations_tuple] = base_distribution
            # Compute difference.
            combination = base_populations_tuple + alternative_populations_tuple
            if combination in min_epsilons_cache:
                min_epsilon = min_epsilons_cache[combination]
            else:
                min_epsilon = crn.CompareDistributions(base_distribution, alternative_distribution, prune=1e-3)
                min_epsilons_cache[combination] = min_epsilon
            # Store difference.
            max_min_epsilon = max(max_min_epsilon, min_epsilon)
        min_epsilons[unsorted_base_populations_tuple] = max_min_epsilon
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
    if args.nspecies == 2:
        values = np.empty(args.nrobots)
        for k, v in min_epsilons.iteritems():
            values[k[1] - 1] = v
        plt.figure()
        plt.plot(np.arange(1, args.nrobots + 1), values, linewidth=2, color='b', label='\epsilon')
        plt.title('N_1 = %d' % args.nrobots)
        plt.xlabel('N_2 (population of the second species)')
        plt.ylabel('Minimum acheivable \epsilon')
        plt.show()
    elif args.nspecies == 3:
        plt.figure()
        x = []
        y = []
        z = []
        for k, v in min_epsilons.iteritems():
            if v == float('inf'):
                continue
            x.append(k[1])
            y.append(k[2])
            z.append(v)
        datapoints = np.array((x, y)).T
        z = np.array(z)
        XI, YI = np.meshgrid(np.arange(1, args.nrobots + 1), np.arange(1, args.nrobots + 1))
        ZI = interpolate.griddata(datapoints, z, (XI, YI), method='linear')
        plt.imshow(ZI, cmap=plt.cm.coolwarm, interpolation='nearest', origin='lower',
                   extent=[0.5, args.nrobots + 0.5, 0.5, args.nrobots + 0.5])
        plt.colorbar()
        plt.title('N_1 = %d' % args.nrobots)
        plt.xlabel('N_2 (population of the second species)')
        plt.ylabel('N_3 (population of the thrid species)')
        plt.show()
    else:
        print 'Cannot plot if --nspecies > 3'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Outputs minimum epsilon for differential privacy')
    parser.add_argument('--stochkit_path', metavar='PATH', action='store', default='./StochKit2.0.11/ssa', help='Path where the stochkit binary is stored')
    parser.add_argument('--alpha', metavar='RATE', type=float, action='store', default=_DEFAULT_ALPHA, help='Alpha parameter (for all species)')
    parser.add_argument('--beta', metavar='RATE', type=float, action='store', default=_DEFAULT_BETA, help='Beta parameter (for all species)')
    parser.add_argument('--duration', metavar='SECONDS', type=float, action='store', default=_DEFAULT_DURATION, help='Run duration')
    parser.add_argument('--nspecies', metavar='N', type=int, action='store', required=True, help='Number of species')
    parser.add_argument('--nrobots', metavar='N', type=int, action='store', required=True, help='Maximum number of robots')
    parser.add_argument('--nruns', metavar='N', type=int, action='store', default=1000, help='Number of StochKit runs')
    parser.add_argument('--cme', action='store_true', help='If set, uses CMEPy instead of StochKit')
    parser.add_argument('--save_epsilons', metavar='FILE', type=str, action='store', help='If set, saves the min-epsilons in the specified file')
    parser.add_argument('--load_epsilons', metavar='FILE', type=str, action='store', help='If set, load the min-epsilons from the specified file')
    parser.add_argument('--plot_first_run', action='store_true', help='If set, plots the first run')

    # Run.
    run(parser.parse_args())
