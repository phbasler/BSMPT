from sympy import *
from sympy.physics.quantum import *
from sympy.printing.cxxcode import cxxcode

import BSMPTTools as BSMPTTools



#parameters
msq, la = symbols('msq lambda', real=True)
params=[msq,la]

#CT params
dmsq, dla = symbols('dmsq dlambda', real=True)
dparams=[dmsq,dla]

#VEVs
v = symbols('v', real=True)

#Higgsfields
rho,eta,zeta,psi = symbols('rho eta zeta psi', real=True)
Higgsfields=[rho,eta,zeta,psi]
CTTadpoles = symbols('dT1:{}'.format(len(Higgsfields)+1))

#doublets
phi = Matrix([[rho+I*eta], [zeta+I*psi]])

#replacements
zeroTempVEV = [(rho,0),(eta,0),(zeta,v),(psi,0)]
fieldsZero = [(x,0) for x in Higgsfields]

phiSq = simplify((Dagger(phi)*phi)[0])
VHiggs = msq/2 * phiSq + la/factorial(4) * phiSq**2


toyModel = BSMPTTools.ModelGenerator(params,dparams,CTTadpoles,Higgsfields,VHiggs)
#toyModel.printTensors()

# get VCT
VCT=toyModel.getVCT()
print("VCT = " + str(VCT))


toyModel.printCTForCPP(zeroTempVEV)