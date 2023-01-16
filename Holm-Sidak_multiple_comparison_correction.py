# Holm-Sidak correction for multiple hypotheses tests
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

# Enter p values data (example below)
p_values = [0.5685, 0.4442, 0.3105, 0.0948, 0.0343]
adjusted_p = []

# Calculate the Holm-Sidak corrected p value
sorted_p_values = sorted(p_values)
for m in range(len(sorted_p_values)):
    if m == 0:
        p_adjust = 1 - (1 - sorted_p_values[m]) ** (len(sorted_p_values))
        adjusted_p.append(p_adjust)
    else:
        p_adjust = max(adjusted_p[m - 1],
                       1 - (1 - sorted_p_values[m]) ** (
                                   len(sorted_p_values) - m + 1))
        adjusted_p.append(p_adjust)

for s in range(len(p_values)):
    print(
        f"The uncorrected: Holm-Sidak corrected p values for #{s + 1} are: "
        f"{sorted_p_values[s] :.5f} : {adjusted_p[s]:.5f}")
