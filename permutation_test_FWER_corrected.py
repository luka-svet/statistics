# Permutation test with correction for multiple hypotheses tests
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


import random

from pandas import *


def draw_samples(pop_un1, pop_tr1, pop_un2, pop_tr2):
    """Draw 4 different populations"""

    # Initialize
    pop_un1_new = []
    pop_tr1_new = []
    pop_un2_new = []
    pop_tr2_new = []

    # Number of observations
    n_pop_un1 = len(pop_un1)
    n_pop_tr1 = len(pop_tr1)

    # Combine untreated populations and treated populations
    pop_comp_un = pop_un1 + pop_un2
    pop_comp_tr = pop_tr1 + pop_tr2

    # Shuffle the combined populations
    random.shuffle(pop_comp_un)
    random.shuffle(pop_comp_tr)

    # Make new random untreated populations and treated populations
    for i in range(n_pop_un1):
        pop_un1_new.append(pop_comp_un[i])

    for j in range(n_pop_tr1):
        pop_tr1_new.append(pop_comp_tr[j])

    for k in range(n_pop_un1, len(pop_comp_un)):
        pop_un2_new.append(pop_comp_un[k])

    for l in range(n_pop_tr1, len(pop_comp_tr)):
        pop_tr2_new.append(pop_comp_tr[l])

    return pop_un1_new, pop_tr1_new, pop_un2_new, pop_tr2_new


def average_population_new(pop_un1_new, pop_tr1_new, pop_un2_new, pop_tr2_new):
    """"Calculate new average population sizes."""

    # Initialize
    avg_pop_un1_new = 0
    avg_pop_tr1_new = 0
    avg_pop_un2_new = 0
    avg_pop_tr2_new = 0

    # Calculate averages
    for observation1 in pop_un1_new:
        avg_pop_un1_new += observation1

    for observation2 in pop_tr1_new:
        avg_pop_tr1_new += observation2

    for observation3 in pop_un2_new:
        avg_pop_un2_new += observation3

    for observation4 in pop_tr2_new:
        avg_pop_tr2_new += observation4

    avg_pop_un1_new_p = avg_pop_un1_new / len(pop_un1_new)
    avg_pop_tr1_new_p = avg_pop_tr1_new / len(pop_tr1_new)
    avg_pop_un2_new_p = avg_pop_un2_new / len(pop_un2_new)
    avg_pop_tr2_new_p = avg_pop_tr2_new / len(pop_tr2_new)

    return avg_pop_un1_new_p, avg_pop_tr1_new_p, avg_pop_un2_new_p, \
        avg_pop_tr2_new_p


def tolerances_new(avg_pop_un1_new_p, avg_pop_tr1_new_p, avg_pop_un2_new_p,
                   avg_pop_tr2_new_p):
    """"Calculate new tolerances."""
    tol1 = avg_pop_tr1_new_p / avg_pop_un1_new_p
    tol2 = avg_pop_tr2_new_p / avg_pop_un2_new_p
    return tol1, tol2


def tolerances_old(pop_un1_old, pop_tr1_old, pop_un2_old, pop_tr2_old):
    """"Calculate old (observed) tolerances."""
    avg_pop_obs = average_population_new(pop_un1_old, pop_tr1_old, pop_un2_old,
                                         pop_tr2_old)
    tolerances_obs = tolerances_new(avg_pop_obs[0], avg_pop_obs[1],
                                    avg_pop_obs[2], avg_pop_obs[3])
    return tolerances_obs


# Main programme

# Data
# ST pop 1 (example below)
# pop_un1_obs = [1000, 2000, 3000]

# pop_tr1_obs = [100, 200, 300]

# ST pop2, 3,...  (example below - enter as many distinct pop. as you have)
# null hypothesis tolerances equal, alt. hypothesis tolerances not equal
# corresponding untreated and treated data has to be in the right order
# pop_un2_obs = [[1000, 2000], [2000, 3000]]
# pop_tr2_obs = [[500, 600], [700, 800]]

start_time = time.time()
# Or import using the csv file
data = read_csv("C:/Users/DrPai/Desktop/data.csv", sep=";")

# converting column data to list,
pop_tr1_obs = data['NUM1'].dropna().tolist()
pop_un1_obs = data['DENOM1'].dropna().tolist()

# import the number of pop_un2_obs and pop_tr2_obs as needed (follow the
# above example)
pop_tr2_obs = [data['NUM2'].dropna().tolist()]
pop_un2_obs = [data['DENOM2'].dropna().tolist()]

# Other parameters
simulations = 30_000

# Number of comparisons vs WT to be done
num_samples = max(len(pop_tr2_obs), len(pop_un2_obs))

# Various lists
tol_old = []
t_distr_all = []
max_t_distr_all = []
t_distr_observed = []
adjusted_p = []
freq_temp = []
frequencies = [0] * num_samples
frequencies_c = [0] * num_samples

for w in range(num_samples):
    # Append the observed tolerances to list tol_old, avoid duplicating WT
    # tolerance
    if not tol_old:  # Check whether the list is empty
        tol_old = list(tolerances_old(pop_un1_obs, pop_tr1_obs, pop_un2_obs[w],
                                      pop_tr2_obs[w]))
    else:  # Append only non-WT tolerance
        tol_old.append(tolerances_old(pop_un1_obs, pop_tr1_obs, pop_un2_obs[w],
                                      pop_tr2_obs[w])[1])
    # Append the observed statistic (for minP-maxT method)
    t_distr_observed.append(abs(tol_old[0] - tol_old[w + 1]))

for x in range(simulations):
    for y in range(num_samples):
        # Draws new (random) samples from the range of observed "populations"
        samples_new = draw_samples(pop_un1_obs, pop_tr1_obs, pop_un2_obs[y],
                                   pop_tr2_obs[y])

        # Calculates the average of new (random) samples from the range of
        # observed "populations"
        avg_pop_new = average_population_new(samples_new[0], samples_new[1],
                                             samples_new[2], samples_new[3])

        # Calculates the tolerance of new (random) samples from the range of
        # observed "populations"
        tol_new = tolerances_new(avg_pop_new[0], avg_pop_new[1],
                                 avg_pop_new[2], avg_pop_new[3])

        # Calculate the two-sided p value
        if (abs(tol_new[0] - tol_new[1])) >= (
                abs(tol_old[0] - tol_old[y + 1])):
            frequencies[y] += 1

        # Calculate the statistic for all the samples from the range of
        # "populations"
        t_distr_all.append(abs(tol_new[0] - tol_new[1]))

    # Save only the maximum (most extreme) value and empty the
    # t_distr_all list
    max_t_distr_all.append(max(t_distr_all))
    t_distr_all.clear()

for q in range(len(frequencies)):
    # Calculate the minP/maxT corrected two-sided p value
    for r in range(len(max_t_distr_all)):
        if max_t_distr_all[r] >= t_distr_observed[q]:
            frequencies_c[q] += 1

# Calculate the Holm-Sidak corrected two-sided p value
sorted_freq = sorted(frequencies)
for m in range(len(sorted_freq)):
    # noinspection PyTypeChecker
    sorted_freq[m] = sorted_freq[m] / simulations
    if m == 0:
        p_adjust = 1 - (1 - sorted_freq[m]) ** (len(sorted_freq))
        adjusted_p.append(p_adjust)
    else:
        p_adjust = max(adjusted_p[m - 1],
                       1 - (1 - sorted_freq[m]) ** (len(sorted_freq) - m + 1))
        adjusted_p.append(p_adjust)

# Find the corresponding two-sided p value for Holm-Sidak and print the results
for g in frequencies:
    freq_temp.append(g / simulations)

for s in range(len(freq_temp)):
    # noinspection PyTypeChecker
    idx = sorted_freq.index(freq_temp[s])
    print(freq_temp[s])
    print(idx)
    print(
        f"The uncorrected: minP/maxT: Holm-Sidak corrected p values for #{s + 1} are: "
        f"{frequencies[s] / simulations:.5f} : "
        f"{frequencies_c[s] / simulations:.5f} : {adjusted_p[idx]:.5f}")
