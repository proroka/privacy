import matplotlib.pyplot as plt
import numpy as np

# Utility class that takes the results from a CRN::Run() and
# shows different types of plots.
class Plotter(object):
    def __init__(self, output_from_run, exclude_list=()):
        self.species, self.timestamps, self.trajectories = output_from_run
        self.nspecies = len(self.species)
        self.exclude_list = exclude_list

    def AverageTrajectory(self):
        mean_trajectory = np.mean(self.trajectories, axis=0)  # First axis is runs.
        std_trajectory = np.std(self.trajectories, axis=0)
        colors = _GetColors(self.nspecies)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in xrange(self.nspecies):
            if self.species[i] in self.exclude_list:
                continue
            plt.fill_between(self.timestamps,
                             mean_trajectory[:, i] + std_trajectory[:, i],
                             mean_trajectory[:, i] - std_trajectory[:, i],
                             facecolor=colors[i], alpha=0.3)
        for i in xrange(self.nspecies):
            if self.species[i] in self.exclude_list:
                continue
            plt.plot(self.timestamps, mean_trajectory[:, i], linewidth=2, color=colors[i], label=self.species[i])
        plt.legend(loc='upper right', shadow=False, fontsize='x-large')
        plt.xlabel('Time [s]')
        plt.ylabel('Population')
        return fig, ax

    def Distributions(self, from_timestamp=0.):
        colors = _GetColors(self.nspecies)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        indices = self.timestamps >= from_timestamp
        min_population = np.max(self.trajectories) + 1
        max_population = -1.
        for i in xrange(self.nspecies):
            if self.species[i] in self.exclude_list:
                continue
            min_population = min(min_population, np.min(self.trajectories[:, indices, i]))
            max_population = max(max_population, np.max(self.trajectories[:, indices, i]))
        for i in xrange(self.nspecies):
            if self.species[i] in self.exclude_list:
                continue
            populations = self.trajectories[:, indices, i].flatten().astype(int)
            counts = np.bincount(populations)
            plt.hist(populations, normed=True, bins=np.arange(min_population - 0.5, max_population + 0.5, 1.),
                     label=self.species[i], color=colors[i], alpha=0.3)
        plt.legend(loc='upper right', shadow=False, fontsize='x-large')
        plt.xlabel('Population')
        plt.ylabel('Probability')
        return fig, ax


def _GetColors(n):
    cm = plt.get_cmap('gist_rainbow')
    return [cm(float(i) / float(n)) for i in range(n)]
