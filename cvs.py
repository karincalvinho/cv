from scipy import stats
import csv
import math
import matplotlib.pyplot as plt
import os
import sys
import numpy
import yaml
import getopt

# process input file
def process_file(input_file, sample_id):
    # extract filename and path from input_files
    path, filename = os.path.split(os.path.abspath(input_file))

    # read input file
    input_data = numpy.genfromtxt(input_file, delimiter=',', skip_header=cfg['skip_header'])
    input_data = numpy.vstack((input_data, input_data[0]))
    v = input_data[:,:1]
    current = input_data[:,1:2]

    # convert Potential vs Hg/Hg0 into Potential vs NHE
    v_nhe = v + cfg['ref'] + .0592 * cfg['pH']

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
    v = v[:endpoint+1, :]
    v_nhe = v_nhe[:endpoint+1, :]
    forward_current = current[:endpoint+1]
    reverse_current = numpy.flipud(current[endpoint:])
    del current

    # find out which vector is the shortest one and trim vectors to that size
    length = min(len(v[:,0]), len(v_nhe[:,0]), len(forward_current[:,0]), len(reverse_current[:,0]))
    v = v[:length]
    v_nhe = v_nhe[:length]
    forward_current = forward_current[:length]
    reverse_current = reverse_current[:length]

    # calculate current density
    current_density = (forward_current + reverse_current) * 1000 / (2 * .196)

    # calculate logarithm of current density
    log10_current_density = numpy.log10(2 * numpy.fabs(current_density))

    # calculate overpotential
    overpotential = v_nhe - 1.23

    # build x/y vectors for linear regression (extracts first quadrant
    # from log10_current_density vs overpotential)
    mask = numpy.logical_and(log10_current_density > -0.5, overpotential > 0.32)
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
        print '%s,error' % input_file
        sys.exit()

    tafel_slope = 1000 * slope

    onset_potential = None
    for (current_density_value, v_nhe_value) in zip(current_density[:,0], v_nhe[:,0]):
        if math.fabs(current_density_value) > math.fabs(cfg['onset']):
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
        'log j (log10 mA/cm2)',
        'Overpotential (V)',
    ]
    output_file = path + '/linear-' + filename
    numpy.savetxt(output_file, linear_part, delimiter=',', header=','.join(column_names),
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

    print '%s,%.3f,%.3f,%.3f,%.3e,%.3f,%.3f' % data_summary



    ################################################################################
    # Print output file
    #

    column_names = [
        'V vs Hg/HgO (V)',
        'i (A)',
        'V vs NHE (V)',
        'j (mA/cm^2)',
        'log j (decade)',
        'Overpotential (V)',
    ]
    output_file = path + '/out-' + filename
    output_data = numpy.hstack((
        v,
        forward_current,
        v_nhe,
        current_density,
        log10_current_density,
        overpotential,
    ))
    numpy.savetxt(output_file, output_data, delimiter=',', header=','.join(column_names),
        comments='')



    ################################################################################
    # Generating Tafel plot
    #

    # chart options
    plt.figure()
    plt.rc('lines', linewidth=2, markersize=4)
    plt.grid(False)
    plt.xlabel('log j (decade)')
    plt.ylabel('Overpotential (V)')

    # plot linear regression
    line_y = []
    for val in x:
        line_y.append(slope * val + intercept)
    plt.plot(x, line_y, '-', color='blue')

    # highlight part of the line between best_i and best_j
    bestline_y = []
    for val in x[best_i:best_j]:
        bestline_y.append(slope * val + intercept)
    plt.plot(x[best_i:best_j], bestline_y, '-', color='red')

    # plot actual data
    plt.plot(x, y, '.', color='black')

    # save Tafel plot to a file
    plt.savefig(path + '/tafel-' + filename + '.png')

    # create multiple line tafel plot
    plt.figure('tafel_all')
    plt.plot(x, y, '.', color='black')
    ax.plot(x, line_y, '-', label=sample_id)
    #plt.plot(x[best_i:best_j], bestline_y, '-', color='#4D4D4D')

    ################################################################################
    # Generating polarization curve
    #

    # chart options
    plt.figure()
    plt.rc('lines', linewidth=2, markersize=4)
    plt.grid(False)
    plt.xlabel('V vs NHE (V)')
    plt.ylabel('Current density (mA/cm2)')
    plt.gca().invert_yaxis()

    # actually plot the data
    plt.plot(v_nhe, current_density)

    # save plot to a file
    plt.savefig(path + '/pc-' + filename + '.png')

def usage():
    print 'Usage: %s <input file>:<sample_id>' % sys.argv[0]
    sys.exit()

# reading script's input from command-line
if len(sys.argv) < 2:
    usage()

config_filename = 'default.yaml'

opts, args = getopt.getopt(sys.argv[1:], 'hc:')
for opt, arg in opts:
    if opt == '-h':
        usage()
    elif opt == '-c':
        config_filename = arg

# reads config file
config_file = open(config_filename, 'r')
cfg = yaml.load(config_file)

# chart options
plt.figure('tafel_all')
ax = plt.gca()
ax.set_color_cycle([
#4D4D4D #(gray)
'#5DA5DA', #(blue)
'#FAA43A', #(orange)
'#60BD68', #(green)
'#F17CB0', #(pink)
'#B2912F', #(brown)
'#B276B2', #(purple)
#DECF3F #(yellow)
'#F15854', #(red)
])
plt.rc('lines', linewidth=2, markersize=4)
plt.grid(False)
plt.xlabel('log j (decade)')
plt.ylabel('Overpotential (V)')

# treat IV data
for argument in args:
    (input_file, sample_id) = argument.split(':')
    process_file(input_file, sample_id)

# save Tafel plot to a file
plt.figure('tafel_all')
plt.legend(loc='lower right')
plt.savefig('tafel_all.png')
