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

import numpy as np
import matplotlib.pyplot as plt
from pandas import *

# Enter data (example below)
# population_d1 = [100, 200, 300]

# population_d2 = [1000, 2000, 3000]

# Or import using the csv file
data = read_csv("C:/data.csv", sep=";")

# converting column data to list
population_d1 = data['NUM'].dropna().tolist()
population_d2 = data['DENOM'].dropna().tolist()


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


# Draw 30000 bootstrap replicates
bs_replicates_array = draw_bs_replicates(population_d1, population_d2, np.mean,
                                         30_000)

# Print empirical mean
print(
    "Empirical mean: " + str(np.mean(population_d1) / np.mean(population_d2)))

# Print the mean of bootstrap replicates
print("Bootstrap replicates mean: " + str(np.mean(bs_replicates_array)))

# Get the corresponding values of 2.5th and 97.5th percentiles
conf_interval = np.percentile(bs_replicates_array, [2.5, 50, 97.5])

# Print the interval
print("The 2.5th, 50th and 97.5th percentile: ", conf_interval)

# Plot the PDF for bootstrap replicates as histogram
plt.hist(bs_replicates_array, bins=200,
         weights=np.ones(len(bs_replicates_array)) / len(bs_replicates_array))

# Showing the related percentiles
plt.axvline(x=np.percentile(bs_replicates_array, [2.5]), ymin=0, ymax=1,
            label='2.5th percentile', c='y')
plt.axvline(x=np.percentile(bs_replicates_array, [50]), ymin=0, ymax=1,
            label='50th percentile', c='g')
plt.axvline(x=np.percentile(bs_replicates_array, [97.5]), ymin=0, ymax=1,
            label='97.5th percentile', c='r')

plt.xlabel("Ratios of means")
plt.ylabel("PDF")
plt.title("Probability Density Function")
plt.legend()
plt.show()
