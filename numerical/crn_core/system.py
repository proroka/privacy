import cmepy.recorder
import numpy as np

from crn import CRN

# Default rate for all species for EXPLORE -> WAIT.
_DEFAULT_ALPHA = 0.2

# Default rate for all species for WAIT -> EXPLORE (without collaboration).
_DEFAULT_BETA = 0.1


def TaskCRNSetup(nrobots, stochkit_binary, alpha=_DEFAULT_ALPHA, beta=_DEFAULT_BETA):
    crn = CRN(stochkit_binary=stochkit_binary)
    nspecies = len(nrobots)
    N = [str(i) for i in xrange(1, nspecies + 1)]  # Species names (i.e., 1, 2, 3, ...).
    # Prepare species alpha/beta.
    if not isinstance(alpha, list) and not isinstance(alpha, tuple):
        alpha = [alpha] * nspecies
    if not isinstance(beta, list) and not isinstance(beta, tuple):
        beta = [beta] * nspecies
    assert len(alpha) == nspecies, 'Size of alpha must be equal to the number of species.'
    assert len(beta) == nspecies, 'Size of beta must be equal to the number of species.'
    # Species go from EXPLORE to WAIT.
    # We will use 'e' for explore and 'w' for wait.
    # A condition like '124w' means that one robot of species 1, 2 and 4 are waiting.
    # TODO: Use a configured alpha and beta for each species.
    waiting_conditions = set()
    for i in xrange(nspecies):
        crn.AddReaction([N[i] + 'e'], [N[i] + 'w'], alpha[i])
        crn.AddReaction([N[i] + 'w'], [N[i] + 'e'], beta[i])
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
                    crn.AddReaction([condition_string, N[i] + 'e'],
                                    [(N[j] + 'e') for j in xrange(nspecies)],
                                    alpha[i])
                else:
                    new_condition_string = ''.join([N[j] for j in new_condition]) + 'w'
                    # Add forward and backward reactions
                    crn.AddReaction([condition_string, N[i] + 'e'], [new_condition_string], alpha[i])
                    crn.AddReaction([new_condition_string], [condition_string, N[i] + 'e'], beta[i])
                    new_waiting_conditions.add(new_condition)
        waiting_conditions = new_waiting_conditions
    # Set population of robots (starting in exploration).
    for i, n in enumerate(nrobots):
        crn.SetInitialPopulation(N[i] + 'e', n)
    return crn


# This function is very hacky. In particular it assumes that the conditions are
# named like: 'Xe', 'Xw', 'XYw', ... and thus uses the length and the last letter
# of the condition to aggregate the data.
def TaskCRNExtractObservableData(output_from_run):
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
