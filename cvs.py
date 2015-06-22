from scipy import stats
import csv
import math
import matplotlib.pyplot as plt
import os
import sys

# reading script's input from command-line
if len(sys.argv) < 2:
    print "Usage: %s <input file>" % sys.argv[0]
    sys.exit()

# input/output filenames
input_file = sys.argv[1]
path, filename = os.path.split(input_file)
output_file = path + '/out-' + filename
chart_file = path + '/chart-' + filename + '.png'

# create data matrix
data = []

# read input file
print ">> Reading input from %s" % input_file
with open(input_file, 'rb') as csvfile:
    rows = csv.reader(csvfile)
    i = 0
    for row in rows:
        i = i + 1
        if i <= 16:
            continue
        row_float = []
        for x in row:
            row_float.append(float(x))
        data.append(row_float)

# convert Potential vs Hg/Hg0 into Potential vs NHE
i = 0
for row in data:
    v = row[0]
    v_nhe = v + 0.140 + .0592*14
    data[i].append(v_nhe)
    i = i + 1

# find inflection point for Potential vs Hg/Hg0
print ">> Finding inflection point"
inflection = None
v_prev = None
i = 0
for (v, _, _) in data:
    if v < v_prev:
        inflection = i - 1
        break
    v_prev = v
    i = i + 1
assert inflection != None
# FIXME: take care of the case that there are multiple cycles in the data
print "   inflection point value is %f" % data[inflection][0]

# calculate current density
i = 0
for row in data[:inflection+1]:
    current = row[1]
    reverse_current = data[-i-1][1]
    data[i].append((current + reverse_current) * 500/ .196)
    i = i + 1

# calculate overpotential
i = 0
for row in data[:inflection+1]:
    v_nhe = row[2]
    o = v_nhe - 1.23
    data[i].append(o)
    i = i + 1

# calculate logarithm of current density
i = 0
for row in data[:inflection+1]:
    current_density = row[3]
    log_cd = math.log10(2 * math.fabs(current_density))
    data[i].append(log_cd)
    i = i + 1

# build x/y vectors for linear regression
x = []
y = []
for row in data[:inflection+1]:
    o = row[4]
    log_cd = row[5]
    if o < 0 or log_cd < 0:
        continue
    x.append(log_cd)
    y.append(o)

# find best linear regression of the curve
print ">> Finding best linear regression of the curve"
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
print "   log j between [%f, %f]" % (x[best_i], x[best_j-1])
print "   Slope: %f" % slope
print "   Tafel slope (mV/dec): %f" % (1000 * slope)
print "   Intercept: %f" % intercept
print "   R-value: %f" % r_value
print "   Standard error: %f" % std_err

# calculate onset potential
onset_potential = None
for row in data[:inflection+1]:
    current_density = row[3]
    if current_density < -.1:
        v_nhe = row[2]
        onset_potential = v_nhe
        break

print ">> Electrochemical activity"
print "   Overpotential at 10mA/cm^2 (V): %f" % (slope + intercept)
print "   Potential at 10mA/cm^2 (V): %f" % (slope + intercept + 1.23)
print "   Onset potential (V): %f" % onset_potential
print "   Exchange current density (mA/cm2): %.3e" % math.pow(10,(-intercept/slope))

# append columns for linear part of tafel plot
for i in range(best_i, best_j):
    data[i].append(x[i])
    data[i].append(y[i])

# write data into output file
print ">> Writing output to %s" % output_file
with open(output_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        'V vs Hg/HgO (V)',
        'i (A)',
        'V vs NHE (V)',
        'j (mA/cm^2)',
        'Overpotential (V)',
        'log j (decade)',
        'Linear log j (decade)',
        'Linear Overpotential (V)',
        ])
    for row in data[:inflection+1]:
        writer.writerow(row)

# generate charts
plt.rc('lines', linewidth=2, markersize=4)

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

# chart options
plt.grid(True)
plt.xlabel('log j (decade)')
plt.ylabel('Overpotential (V)')

# save chart to a file
print ">> Saving chart to %s" % chart_file
plt.savefig(chart_file)
