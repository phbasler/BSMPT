from sympy import symbols, Matrix, simplify, factorial, I, sqrt
from sympy.physics.quantum import Dagger


import argparse as argparse
import ModelGenerator as ModelGenerator


# SM paramters

Cg = symbols("C_g", real=True)
Cgs = symbols("C_gs", real=True)
sigma0 = Matrix([[1, 0], [0, 1]])
sigma1 = Matrix([[0, 1], [1, 0]])
sigma2 = Matrix([[0, -I], [I, 0]])
sigma3 = Matrix([[1, 0], [0, -1]])
m_electron = symbols("C_MassElectron", real=True)
m_mu = symbols("C_MassMu", real=True)
m_tau = symbols("C_MassTau", real=True)
m_up = symbols("C_MassUp", real=True)
m_charm = symbols("C_MassCharm", real=True)
m_top = symbols("C_MassTop", real=True)
m_down = symbols("C_MassDown", real=True)
m_strange = symbols("C_MassStrange", real=True)
m_bottom = symbols("C_MassBottom", real=True)

# CKM Matrix
Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb = symbols(
    "Vud Vus Vub Vcd Vcs Vcb Vtd Vts Vtb"
)
VCKM = Matrix([[Vud, Vus, Vub], [Vcd, Vcs, Vcb], [Vtd, Vts, Vtb]])

# parameters
msq, la = symbols("msq lambda", real=True)
params = [msq, la]

# CT params
dmsq, dla = symbols("dmsq dlambda", real=True)
dparams = [dmsq, dla]

# VEVs
v = symbols("v", real=True)

# Higgsfields
rho, eta, zeta, psi = symbols("rho eta zeta psi", real=True)
Higgsfields = [rho, eta, zeta, psi]
CTTadpoles = symbols("dT1:{}".format(len(Higgsfields) + 1), real=True)

# doublets
phi = Matrix([[rho + I * eta], [zeta + I * psi]]) * 1 / sqrt(2)

# replacements
zeroTempVEV = [(rho, 0), (eta, 0), (zeta, v), (psi, 0)]
finiteTempVEV = zeroTempVEV
fieldsZero = [(x, 0) for x in Higgsfields]

phiSq = simplify((Dagger(phi) * phi)[0])
VHiggs = msq / 2 * phiSq + la / factorial(4) * phiSq**2

# Set Gauge fields
W1, W2, W3, B0 = symbols("W1 W2 W3 B0", real=True)


Dmu = (
    -I * Cg / 2 * (sigma1 * W1 + sigma2 * W2 + sigma3 * W3) - I * Cgs / 2 * sigma0 * B0
)
VGauge = simplify(Dagger(Dmu * phi) * (Dmu * phi))[0, 0]

# Generate Lepton Potentials
NuL = symbols("veL vmuL vtauL", real=True)
ER = symbols("eR muR tauR", real=True)
EL = symbols("eL muL tauL", real=True)

ye = sqrt(2) * m_electron / v
ymu = sqrt(2) * m_mu / v
ytau = sqrt(2) * m_tau / v

PiLep = Matrix([[ye, 0, 0], [0, ymu, 0], [0, 0, ytau]])

VFLep = 0
for i in range(len(NuL)):
    for j in range(len(ER)):
        VFLep += (NuL[i] * PiLep[i, j] * ER[j]) * phi[0]

for i in range(len(EL)):
    for j in range(len(ER)):
        VFLep += (EL[i] * PiLep[i, j] * ER[j]) * phi[1]

VFLep = simplify(VFLep)
LepBase = symbols('eL eR muL muR tauL tauR veL vmuL vtauL', real=True)

# Generate Quark Potentials
UL = symbols("uL cL tL", real=True)
DL = symbols("dL sL bL", real=True)
UR = symbols("uR cR tR", real=True)
DR = symbols("dR sR bR", real=True)
QuarkBase = UL + DL + UR + DR

yb = sqrt(2) * m_bottom / v
yc = sqrt(2) * m_charm / v
yd = sqrt(2) * m_down / v
ys = sqrt(2) * m_strange / v
yt = sqrt(2) * m_top / v
yu = sqrt(2) * m_up / v

DownCoupling = Matrix([[yd, 0, 0], [0, ys, 0], [0, 0, yb]])
UpCoupling = Matrix([[yu, 0, 0], [0, yc, 0], [0, 0, yt]])

ULVector = Matrix([[x] for x in UL])
DLVector = Matrix([[x] for x in DL])
URVector = Matrix([[x] for x in UR])
DRVector = Matrix([[x] for x in DR])

VQuark = ULVector.transpose() * VCKM * DownCoupling * DRVector * phi[0]
VQuark += DLVector.transpose() * DownCoupling * DRVector * phi[0]
VQuark += ULVector.transpose() * UpCoupling * URVector * phi[1].conjugate()
VQuark += (
    -DLVector.transpose() * Dagger(VCKM) * UpCoupling * URVector * phi[1].conjugate()
)
VQuark = simplify(VQuark[0, 0])

# Generate the model
toyModel = ModelGenerator.ModelGenerator(
    params, dparams, CTTadpoles, Higgsfields, VHiggs, zeroTempVEV, finiteTempVEV
)
toyModel.setGauge([W1, W2, W3, B0], VGauge)
toyModel.setLepton(LepBase, VFLep)
toyModel.setQuark(QuarkBase, VQuark)


parser = argparse.ArgumentParser()
parser.add_argument(
    "-s",
    "--show",
    choices=["ct", "tensor", "treeSimpl", "CTSimpl"],
    required=True,
    help="The part of the model to be printed",
)

if __name__ == "__main__":
    args = parser.parse_args()
    method = args.show

    printCT = method == "ct"
    printTensors = method == "tensor"

    if printCT:
        print("//Begin CT Calculation")
        toyModel.printCTForCPP()
        print("//End CT Calculation")

        print("//Begin CT Order for set_CT_Pot_Par")
        toyModel.printCTOrder()
        print("//End CT Order")

    if printTensors:
        toyModel.printModelToCPP()

    if method == "treeSimpl":
        toyModel.printTreeSimplified()

    if method == "CTSimpl":
        toyModel.printVCTSimplified()
