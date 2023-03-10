# Bootstrap test for ratios with correction for multiple hypotheses tests
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
from pandas import *

# Import using the csv file
data = read_csv("C:/data.csv", sep=";")

# converting column data to list
population_tr = [data['TR1'].dropna().tolist(), data['TR2'].dropna().tolist(),
                 data['TR3'].dropna().tolist(), data['TR4'].dropna().tolist(),
                 data['TR5'].dropna().tolist(), data['TR6'].dropna().tolist(),
                 data['TR7'].dropna().tolist()]
population_un = [data['UN1'].dropna().tolist(), data['UN2'].dropna().tolist(),
                 data['UN3'].dropna().tolist(), data['UN4'].dropna().tolist(),
                 data['UN5'].dropna().tolist(), data['UN6'].dropna().tolist(),
                 data['UN7'].dropna().tolist()]


def draw_bs_replicates(pop_tr1, pop_un1, pop_tr2, pop_un2,
                       size):
    """Creates bootstrap samples, computes replicates and returns
    replicates array."""
    # Create an empty array to store replicates
    bs_replicates = np.empty(size)

    # Create bootstrap replicates as much as size
    for i in range(size):
        # Create bootstrap samples
        bs_sample_tr1 = np.random.choice(pop_tr1, size=len(pop_tr1))
        bs_sample_un1 = np.random.choice(pop_un1, size=len(pop_un1))
        bs_sample_tr2 = np.random.choice(pop_tr2, size=len(pop_tr2))
        bs_sample_un2 = np.random.choice(pop_un2, size=len(pop_un2))
        # Get bootstrap replicate and append to bs_replicates
        bs_replicates[i] = observed_x_statistic(bs_sample_tr1, bs_sample_un1,
                                                bs_sample_tr2, bs_sample_un2)

    return bs_replicates


def observed_x_statistic(tr_data1, un_data1, tr_data2, un_data2):
    """Computes x statistic of the observed data"""
    mean_tr1 = np.mean(tr_data1)
    mean_un1 = np.mean(un_data1)
    mean_tr2 = np.mean(tr_data2)
    mean_un2 = np.mean(un_data2)
    obs_x_statistic = (mean_tr1 / mean_un1) - (mean_tr2 / mean_un2)

    return obs_x_statistic


# General variables
simulations = 30_000

# Number of comparisons vs WT to be done
num_samples = len(population_tr) - 1

# Various lists
observed_x = []
bs_replicates_all = []
frequencies = [
                  1] * num_samples  # Better than starting with 0 (see https://stats.stackexchange.com/questions/92542/how-to-perform-a-bootstrap-test-to-compare-the-means-of-two-samples)
adjusted_p = []
freq_temp = []

for w in range(num_samples):
    # Append the observed statistic
    observed_x.append(observed_x_statistic(population_tr[0], population_un[0],
                                           population_tr[w + 1],
                                           population_un[w + 1]))

    # Draw n bootstrap replicates for all the samples
    bs_replicates_all.append(
        draw_bs_replicates(population_tr[0], population_un[0],
                           population_tr[w + 1], population_un[w + 1],
                           simulations))
    for y in range(simulations):
        if observed_x[w] >= 0:
            if bs_replicates_all[w][y] <= 0:
                frequencies[w] += 1
        elif observed_x[w] < 0:
            if bs_replicates_all[w][y] >= 0:
                frequencies[w] += 1

# Calculate the Holm-Sidak corrected two-sided p value
sorted_freq = sorted(frequencies)
for m in range(len(sorted_freq)):
    # noinspection PyTypeChecker
    sorted_freq[m] = sorted_freq[m] / (
            simulations + 1)  # Better than tarting with 1, see above
    if m == 0:
        p_adjust = 1 - (1 - sorted_freq[m]) ** (len(sorted_freq))
        adjusted_p.append(p_adjust)
    else:
        p_adjust = max(adjusted_p[m - 1],
                       1 - (1 - sorted_freq[m]) ** (len(sorted_freq) - m + 1))
        adjusted_p.append(p_adjust)

# Find the corresponding two-sided p value for Holm-Sidak and print the results
for g in frequencies:
    freq_temp.append(g / (simulations + 1))

for s in range(len(freq_temp)):
    # noinspection PyTypeChecker
    idx = sorted_freq.index(freq_temp[s])
    print(
        f"The uncorrected : Holm-Sidak corrected p values for #{s + 1} are: "
        f"{frequencies[s] / (simulations + 1):.5f} : "
        f"{adjusted_p[idx]:.5f}")
