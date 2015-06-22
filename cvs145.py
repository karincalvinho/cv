from scipy import stats
import csv
import math
import matplotlib.pyplot as plt
import os
import sys
import numpy

# reading script's input from command-line
if len(sys.argv) < 2:
    print "Usage: %s <input file>" % sys.argv[0]
    sys.exit()

# extract filename and path from input_file
input_file = sys.argv[1]
path, filename = os.path.split(input_file)

# read input file
data = numpy.genfromtxt(input_file, delimiter=',', skip_header=16)
data = numpy.vstack((data, data[0]))
v = data[:,:1]
current = data[:,1:2]

# convert Potential vs Hg/Hg0 into Potential vs NHE
v_nhe = v + 0.145 + .0592*14
data = numpy.append(data, v_nhe, 1)

# find cycle endpoint
endpoint = None
v_prev = None
i = 0
for v_now in v[:,0]:
    if v_now < v_prev:
        endpoint = i - 1
        break
    v_prev = v_now
    i = i + 1
assert endpoint != None
# FIXME: take care of the case that there are multiple cycles in the data

# trimming data, excluding rows after endpoint from output
data = data[:endpoint+1, :]

# calculate current density
forward_current = current[:endpoint+1]
reverse_current = numpy.flipud(current[endpoint:])
current_density = (forward_current + reverse_current) * 1000 / (2 * .196)
data = numpy.append(data, current_density, 1)

# calculate logarithm of current density
log10_current_density = numpy.log10(2 * numpy.fabs(current_density))
data = numpy.append(data, log10_current_density, 1)

# calculate overpotential
overpotential = v_nhe[:endpoint+1] - 1.23
data = numpy.append(data, overpotential, 1)

# build x/y vectors for linear regression (extracts first quadrant
# from log10_current_density vs overpotential)
mask = numpy.logical_and(log10_current_density > 0, overpotential > 0)
x = log10_current_density[mask]
y = overpotential[mask]



################################################################################
# Linear regression of the curve
#

# find best linear regression of the curve
best_r_value = 0
best_i = None
best_j = None
for i in range(0, len(x)):
    for j in range(i+1, len(x)+1):
        if x[j-1] - x[i] < 1:
            continue

        slope, intercept, r_value, p_value, std_err = stats.linregress(
                x[i:j], y[i:j])

        if r_value > best_r_value:
            best_r_value = r_value
            best_i = i
            best_j = j

# prints the best linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(
    x[best_i:best_j], y[best_i:best_j])

if best_i == None or best_j == None:
    print "%s,error" % input_file
    sys.exit()

tafel_slope = 1000 * slope

onset_potential = None
for row in data:
    current_density_value = row[3]
    if current_density_value < -1:
        v_nhe_value = row[2]
        onset_potential = v_nhe_value
        break

potential_at_10 = slope + intercept + 1.23

exchange_current_density = math.pow(10,(-intercept/slope))

# write linear part of tafel plot into a separate file
linear_part = numpy.hstack((
  numpy.vstack(x[best_i:best_j]),
  numpy.vstack(y[best_i:best_j])
  ))
column_names = [
    'Overpotential (V)',
    'log j (decade)',
]
output_file = path + '/linear-' + filename
numpy.savetxt(output_file, data, delimiter=',', header=','.join(column_names),
    comments='')



################################################################################
# Print data summary
#

data_summary = (
  input_file,
  tafel_slope,
  onset_potential,
  potential_at_10,
  exchange_current_density,
  intercept,
  r_value,
)

print "%s,%.3f,%.3f,%.3f,%.3e,%.3f,%.3f" % data_summary



################################################################################
# Print output file
#

column_names = [
    'V vs Hg/HgO (V)',
    'i (A)',
    'V vs NHE (V)',
    'j (mA/cm^2)',
    'Overpotential (V)',
    'log j (decade)',
]
output_file = path + '/out-' + filename
numpy.savetxt(output_file, data, delimiter=',', header=','.join(column_names),
    comments='')



################################################################################
# Generating Tafel plot
#

# chart options
plt.clf()
plt.rc('lines', linewidth=2, markersize=4)
plt.grid(True)
plt.xlabel('log j (decade)')
plt.ylabel('Overpotential (V)')

# plot linear regression
line_y = []
for val in x:
    line_y.append(slope * val + intercept)
plt.plot(x, line_y, '-', color="blue")

# highlight part of the line between best_i and best_j
line_y = []
for val in x[best_i:best_j]:
    line_y.append(slope * val + intercept)
plt.plot(x[best_i:best_j], line_y, '-', color="red")

# plot actual data
plt.plot(x, y, '.k')

# save Tafel plot to a file
plt.savefig(path + '/tafel-' + filename + '.png')



################################################################################
# Generating polarization curve
#

# chart options
plt.clf()
plt.rc('lines', linewidth=2, markersize=4)
plt.grid(True)
plt.xlabel('V vs NHE (V)')
plt.ylabel('Current density (mA/cm2)')
plt.gca().invert_yaxis()

# actually plot the data
plt.plot(v_nhe[:701], current_density[:701])

# save plot to a file
plt.savefig(path + '/pc-' + filename + '.png')
