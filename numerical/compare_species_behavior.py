import argparse
import itertools
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import pickle
from scipy import interpolate

import crn_core

# Default rate for all species for EXPLORE -> WAIT.
_DEFAULT_ALPHA = 1.0

# Default rate for all species for WAIT -> EXPLORE (without collaboration).
_DEFAULT_BETA = 1.0

# Default duration
_DEFAULT_DURATION = 5.

# Number of points per axis in the graph.
_NPOINTS = 10

# To only show the first call to RunSystem.
first_run = True


def RunSystem(populations, alphas, betas, args):
    global first_run

    # Ignore collaboration rate.
    crn_model = crn_core.TaskCRNSetup(populations, args.stochkit_path, alpha=alphas, beta=betas)
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
        print 'Showing plot for alpha:', alphas, 'beta:', betas
        p = crn_core.Plotter(aggregated_output)
        p.AverageTrajectory()
        p.Distributions(from_timestamp=args.duration / 2.0)
        plt.show()
        first_run = False
    return crn_core.BuildDistribution(aggregated_output, from_timestamp=args.duration / 2.0)


def GetAllEpsilons(args):
    assert args.sweep_type == 'alpha' or args.sweep_type == 'beta', 'Only types "alpha" and "beta" are supported.'
    # Base population (all species have the same number of robots).
    base_populations = np.array([args.nrobots] * args.nspecies)
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
    if args.sweep_type == 'alpha':
        parameters = np.linspace(0., args.alpha * 2., args.npoints + 1).tolist()[1:]  # Ignore zero.
        default = args.alpha
    else:
        parameters = np.linspace(0., args.beta * 2., args.npoints + 1).tolist()[1:]  # Ignore zero.
        default = args.beta
    parameters = [(default,) + x for x in itertools.combinations_with_replacement(parameters, args.nspecies - 1)]
    # Loop through all alpha parameters.
    min_epsilons = {}
    for parameter in parameters:
        print 'Analyzing', parameter
        if args.sweep_type == 'alpha':
            base_distribution = RunSystem(base_populations, parameter, args.beta, args)
        else:
            base_distribution = RunSystem(base_populations, args.alpha, parameter, args)
        # Alternative database.
        max_min_epsilon = 0.
        for alternative_populations in alternative_population_offsets:
            if args.sweep_type == 'alpha':
                alternative_distribution = RunSystem(alternative_populations + base_populations, parameter, args.beta, args)
            else:
                alternative_distribution = RunSystem(alternative_populations + base_populations, args.alpha, parameter, args)
            # Compute difference.
            min_epsilon = crn_core.CompareDistributions(base_distribution, alternative_distribution, smooth=1e-3)
            # Store difference.
            max_min_epsilon = max(max_min_epsilon, min_epsilon)
        # Due to symmetry all permutations of the parameters (except the first element) should be equivalent.
        for params in ((parameter[0],) + p for p in itertools.permutations(parameter[1:])):
            min_epsilons[params] = max_min_epsilon
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
        values = np.array(sorted([(k[1], v) for k, v in min_epsilons.iteritems()]))
        plt.figure()
        plt.plot(values[:, 0], values[:, 1], linewidth=2, color='b', label='\epsilon')
        plt.title('N_1 = N_2 = %d, %s_1 = %g, %s = %g' % (
            args.nrobots, args.sweep_type, args.alpha if args.sweep_type == 'alpha' else args.beta,
            'beta' if args.sweep_type == 'alpha' else 'alpha', args.beta if args.sweep_type == 'alpha' else args.alpha))
        plt.xlabel('%s_2' % args.sweep_type)
        plt.ylabel('Leakage')
        plt.show()
    elif args.nspecies == 3:
        fig = plt.figure(figsize=(5, 4))
        x = []
        y = []
        z = []
        parameters = set()
        for k, v in min_epsilons.iteritems():
            parameters.add(k[1])
            parameters.add(k[2])
            if v == float('inf'):
                continue
            x.append(k[1])
            y.append(k[2])
            z.append(v)
        datapoints = np.array((x, y)).T
        parameters = np.array(sorted(list(parameters)))
        z = np.array(z)
        XI, YI = np.meshgrid(parameters, parameters)
        ZI = interpolate.griddata(datapoints, z, (XI, YI), method='linear')
        dx = (np.max(parameters) - np.min(parameters)) / (2. * float(len(parameters) - 1))
        if args.sweep_type == 'alpha':
            clim = None # (0., .6)
        else:
            clim = None # (0.04, 0.08)
        norm = colors.PowerNorm(gamma=1)
        cmap = plt.get_cmap('RdPu')
        plt.imshow(ZI, cmap=cmap, interpolation='nearest', origin='lower',
                   clim=clim, norm=norm,
                   extent=[np.min(parameters) - dx, np.max(parameters) + dx, np.min(parameters) - dx, np.max(parameters) + dx])
        plt.colorbar()
        plt.title('N_1 = N_2 = N_3 = %d, %s_1 = %g, %s = %g' % (
            args.nrobots, args.sweep_type, args.alpha if args.sweep_type == 'alpha' else args.beta,
            'beta' if args.sweep_type == 'alpha' else 'alpha', args.beta if args.sweep_type == 'alpha' else args.alpha))
        plt.xlabel('%s_2' % args.sweep_type)
        plt.ylabel('%s_3' % args.sweep_type)
        plt.show()
    else:
        print 'Cannot plot if --nspecies > 3'
    if args.save_plot:
        fig.savefig(args.save_plot)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Outputs minimum epsilon for differential privacy')
    parser.add_argument('--stochkit_path', metavar='PATH', action='store', default='./StochKit2.0.11/ssa', help='Path where the stochkit binary is stored')
    parser.add_argument('--duration', metavar='SECONDS', type=float, action='store', default=_DEFAULT_DURATION, help='Run duration')
    parser.add_argument('--alpha', metavar='RATE', type=float, action='store', default=_DEFAULT_ALPHA, help='Average value of alpha parameter')
    parser.add_argument('--beta', metavar='RATE', type=float, action='store', default=_DEFAULT_BETA, help='Average value of beta parameter')
    parser.add_argument('--nspecies', metavar='N', type=int, action='store', required=True, help='Number of species')
    parser.add_argument('--nrobots', metavar='N', type=int, action='store', required=True, help='Maximum number of robots (for all species)')
    parser.add_argument('--nruns', metavar='N', type=int, action='store', default=1000, help='Number of StochKit runs')
    parser.add_argument('--npoints', metavar='N', type=int, action='store', default=_NPOINTS, help='Number datapoints')
    parser.add_argument('--cme', action='store_true', help='If set, uses CMEPy instead of StochKit')
    parser.add_argument('--save_epsilons', metavar='FILE', type=str, action='store', help='If set, saves the min-epsilons in the specified file')
    parser.add_argument('--load_epsilons', metavar='FILE', type=str, action='store', help='If set, load the min-epsilons from the specified file')
    parser.add_argument('--sweep_type', metavar='{alpha | beta}', type=str, action='store', default='alpha', help='Either "alpha" or "beta"')
    parser.add_argument('--save_plot', metavar='FILE', type=str, action='store', help='If set, saves the plot in the specified file')
    parser.add_argument('--plot_first_run', action='store_true', help='If set, plots the first run')

    # Run.
    run(parser.parse_args())
