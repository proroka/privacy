import argparse
import itertools
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import pickle
from scipy import interpolate

import crn_core

# Default rate for all species for EXPLORE -> WAIT.
_DEFAULT_ALPHA = 0.2

# Default rate for all species for WAIT -> EXPLORE (without collaboration).
_DEFAULT_BETA = 0.1

# Default duration
_DEFAULT_DURATION = 5.

# To only show the first call to RunSystem.
first_run = True

# nrobots going from 1 to _MAX_NROBOTS_FACTOR*nrobots.
_MAX_NROBOTS_FACTOR = 2


def RunSystem(populations, args):
    global first_run

    # Ignore collaboration rate.
    crn_model = crn_core.TaskCRNSetup(populations, args.stochkit_path,  alpha=args.alpha, beta=args.beta)
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
    aggregated_output = crn_core.TaskCRNExtractObservableData(output)
    if args.plot_first_run and first_run:
        print 'Showing plot for the population:', populations
        p = crn_core.Plotter(aggregated_output)
        p.AverageTrajectory()
        p.Distributions(from_timestamp=args.duration / 2.0)
        plt.show()
        first_run = False
    return crn_core.BuildDistribution(aggregated_output, from_timestamp=args.duration / 2.0)


def GetAllEpsilons(args):
    # 1 species is fixed at args.nrobots, while the others vary in numbers.
    populations_to_test = itertools.combinations_with_replacement(xrange(1, args.nrobots * _MAX_NROBOTS_FACTOR + 1), args.nspecies - 1)
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
    # Since the system is symmetric, we store in sorted population numbers.
    # So that we do not test both (N_1, N_2, N_3) and (N_1, N_3, N_2).
    min_epsilons = {}
    for populations in populations_to_test:
        print 'Analyzing population:', populations
        base_distribution = RunSystem(populations, args)
        # Alternative database.
        max_min_epsilon = 0.
        for alternative_populations in alternative_population_offsets:
            alternative_distribution = RunSystem(alternative_populations + populations, args)
            # Compute difference.
            min_epsilon = crn_core.CompareDistributions(base_distribution, alternative_distribution)
            # Store difference.
            max_min_epsilon = max(max_min_epsilon, min_epsilon)
        # Due to symmetry all permutations of the parameters (except the first element) should be equivalent.
        for pop in ((populations[0],) + p for p in itertools.permutations(populations[1:])):
            min_epsilons[pop] = max_min_epsilon
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
        values = np.empty(args.nrobots * _MAX_NROBOTS_FACTOR)
        for k, v in min_epsilons.iteritems():
            values[k[1] - 1] = v
        fig = plt.figure()
        plt.plot(np.arange(1, args.nrobots * _MAX_NROBOTS_FACTOR + 1), values, linewidth=2, color='b', label='\epsilon')
        plt.title('N_1 = %d' % args.nrobots)
        plt.xlabel('N_2 (population of the second species)')
        plt.ylabel('Leakage')
        plt.show()
    elif args.nspecies == 3:
        fig = plt.figure(figsize=(5, 4))
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
        XI, YI = np.meshgrid(np.arange(1, args.nrobots * _MAX_NROBOTS_FACTOR + 1), np.arange(1, args.nrobots * _MAX_NROBOTS_FACTOR + 1))
        ZI = interpolate.griddata(datapoints, z, (XI, YI), method='linear')
        clim = (0.,10.)
        norm = colors.PowerNorm(gamma=1.)
        cmap = plt.get_cmap('RdPu')
        plt.imshow(ZI, cmap=cmap, interpolation='nearest', origin='lower',
                   clim=clim, norm=norm,
                   extent=[0.5, args.nrobots * _MAX_NROBOTS_FACTOR + 0.5, 0.5, args.nrobots * _MAX_NROBOTS_FACTOR + 0.5])
        plt.colorbar()
        plt.title('N_1 = %d' % args.nrobots)
        plt.xlabel('N_2 (population of the second species)')
        plt.ylabel('N_3 (population of the thrid species)')
        plt.show()
    else:
        print 'Cannot plot if --nspecies > 3'
    if args.save_plot:
        fig.savefig(args.save_plot)


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
    parser.add_argument('--save_plot', metavar='FILE', type=str, action='store', help='If set, saves the plot in the specified file')
    parser.add_argument('--plot_first_run', action='store_true', help='If set, plots the first run')

    # Run.
    run(parser.parse_args())
