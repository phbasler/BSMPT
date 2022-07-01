from pyBSMPT import BSMPT

# set input arguments
inputfile = '../../example/C2HDM_Input.dat'
bsmptdir = '../../build'
outputdir = 'BSMPT_output/test.csv'
model = 'c2hdm'

# load point
point = BSMPT(model, inputfile, bsmptdir)

# call BSMPT
firstline = 2
secondline = 2
point.calc_BSMPT(outputdir, firstline, secondline)

# get PT strength
xic = point.get_value(outputdir, 'omega_c/T_c')

print(f'The strength of the EWPT is calculated to be xic = {xic[0]}.')