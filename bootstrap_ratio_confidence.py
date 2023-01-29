# Calculation of confidence intervals for ratios of means via bootstrapping
#
# Copyright 2022-2023 Luka Svet <luka.svet@kuleuven.be>
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt
import numpy as np
from pandas import *

# Or import using the csv file
data = read_csv("C:/data.csv", sep=";")

# converting column data to list
population_d1 = [data['TR1'].dropna().tolist(), data['TR2'].dropna().tolist(),
                 data['TR3'].dropna().tolist(), data['TR4'].dropna().tolist(),
                 data['TR5'].dropna().tolist(), data['TR6'].dropna().tolist(),
                 data['TR7'].dropna().tolist()]
population_d2 = [data['UN1'].dropna().tolist(), data['UN2'].dropna().tolist(),
                 data['UN3'].dropna().tolist(), data['UN4'].dropna().tolist(),
                 data['UN5'].dropna().tolist(), data['UN6'].dropna().tolist(),
                 data['UN7'].dropna().tolist()]


def draw_bs_replicates(data1, data2, func, size):
    """Creates bootstrap samples, computes replicates and returns
    replicates array """
    # Create an empty array to store replicates
    bs_replicates = np.empty(size)

    # Create bootstrap replicates as much as size
    for i in range(size):
        # Create bootstrap samples
        bs_sample1 = np.random.choice(data1, size=len(data1))
        bs_sample2 = np.random.choice(data2, size=len(data2))
        # Get bootstrap replicate and append to bs_replicates
        bs_replicates[i] = func(bs_sample1) / func(bs_sample2)

    return bs_replicates


bs_replicates_array = []
num_samples = len(population_d1)

for w in range(num_samples):
    # Draw 30000 bootstrap replicates for all the samples
    bs_replicates_array.append(
        draw_bs_replicates(population_d1[w], population_d2[w], np.mean,
                           30_000))

    # Print sample number
    print("--------------------" + str(w + 1) + "--------------------")

    # Print the empirical mean
    print(
        "Empirical mean: " + str(
            np.mean(population_d1[w]) / np.mean(population_d2[w])))

    # Print the mean of bootstrap replicates
    print("Bootstrap replicates mean: " + str(np.mean(bs_replicates_array[w])))

    # Get the corresponding values of 2.5th and 97.5th percentiles
    conf_interval = np.percentile(bs_replicates_array[w], [2.5, 50, 97.5])

    # Print the interval
    print("The 2.5th, 50th and 97.5th percentile: ", conf_interval)

    # End block
    print("----------------------------------------\n")

    # Plot the PDF for bootstrap replicates as histogram
    plt.hist(bs_replicates_array[w], bins=200,
             weights=np.ones(len(bs_replicates_array[w])) / len(
                 bs_replicates_array))

    # Showing the related percentiles
    plt.axvline(x=np.percentile(bs_replicates_array[w], [2.5]), ymin=0, ymax=1,
                label='2.5th percentile', c='y')
    plt.axvline(x=np.percentile(bs_replicates_array[w], [50]), ymin=0, ymax=1,
                label='50th percentile', c='g')
    plt.axvline(x=np.percentile(bs_replicates_array[w], [97.5]), ymin=0,
                ymax=1,
                label='97.5th percentile', c='r')

    plt.xlabel("Ratios of means")
    plt.ylabel("PDF")
    plt.title("Probability Density Function")
    plt.legend()
    plt.show()
