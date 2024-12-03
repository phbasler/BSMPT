from sympy import symbols, Matrix, simplify, I, sqrt, im, conjugate, expand
from sympy.physics.quantum import Dagger
from enum import Enum

import ModelGenerator as ModelGenerator
import argparse


class TwoHDMType(Enum):
    TypeI = (0,)
    TypeII = (1,)
    Flipped = (2,)
    LeptonSpecific = 3


Chosen2HDMType = TwoHDMType.TypeI


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
m11sq, m22sq, Rem12sq, Imm12sq = symbols("m11sq m22sq Rem12sq Imm12sq", real=True)
lambda1 = symbols("lambda1", real=True)
lambda2 = symbols("lambda2", real=True)
lambda3 = symbols("lambda3", real=True)
lambda4 = symbols("lambda4", real=True)
Relambda5 = symbols("Relambda5", real=True)
Imlambda5 = symbols("Imlambda5", real=True)
Relambda6 = symbols("Relambda6", real=True)
Imlambda6 = symbols("Imlambda6", real=True)
Relambda7 = symbols("Relambda7", real=True)
Imlambda7 = symbols("Imlambda7", real=True)

m12sq = Rem12sq + I * Imm12sq
lambda5 = Relambda5 + I * Imlambda5
lambda6 = Relambda6 + I * Imlambda6
lambda7 = Relambda7 + I * Imlambda7

params = [
    lambda1,
    lambda2,
    lambda3,
    lambda4,
    Relambda5,
    Imlambda5,
    Relambda6,
    Imlambda6,
    Relambda7,
    Imlambda7,
    m11sq,
    m22sq,
    Rem12sq,
    Imm12sq,
]

# CT params
dm11sq, dm22sq, dRem12sq, dImm12sq = symbols(
    "dm11sq dm22sq dRem12sq dImm12sq", real=True
)
dlambda1 = symbols("dlambda1", real=True)
dlambda2 = symbols("dlambda2", real=True)
dlambda3 = symbols("dlambda3", real=True)
dlambda4 = symbols("dlambda4", real=True)
dRelambda5 = symbols("dRelambda5", real=True)
dImlambda5 = symbols("dImlambda5", real=True)
dRelambda6 = symbols("dRelambda6", real=True)
dImlambda6 = symbols("dImlambda6", real=True)
dRelambda7 = symbols("dRelambda7", real=True)
dImlambda7 = symbols("dImlambda7", real=True)
dparams = [
    dlambda1,
    dlambda2,
    dlambda3,
    dlambda4,
    dRelambda5,
    dImlambda5,
    dRelambda6,
    dImlambda6,
    dRelambda7,
    dImlambda7,
    dm11sq,
    dm22sq,
    dRem12sq,
    dImm12sq,
]

# VEVs
v1 = symbols("v1", real=True)
v2 = symbols("v2", real=True)
w1 = symbols("w1", real=True)
w2 = symbols("w2", real=True)
wCP = symbols("wCP", real=True)
wCB = symbols("wCB", real=True)

# Higgsfields
rho1, eta1, zeta1, psi1 = symbols("rho1 eta1 zeta1 psi1", real=True)
rho2, eta2, zeta2, psi2 = symbols("rho2 eta2 zeta2 psi2", real=True)
Higgsfields = [rho1, eta1, rho2, eta2, zeta1, psi1, zeta2, psi2]
CTTadpoles = symbols("dT1:{}".format(len(Higgsfields) + 1), real=True)

# doublets
phi1 = (
    Matrix(
        [[Higgsfields[0] + I * Higgsfields[1]], [Higgsfields[4] + I * Higgsfields[5]]]
    )
    * 1
    / sqrt(2)
)
phi2 = (
    Matrix(
        [[Higgsfields[2] + I * Higgsfields[3]], [Higgsfields[6] + I * Higgsfields[7]]]
    )
    * 1
    / sqrt(2)
)

# replacements
higgsvevAtZeroTemp = [0, 0, 0, 0, v1, 0, v2, 0]
higgsVEVAtFiniteTemp = [0, 0, wCB, 0, w1, 0, w2, wCP]
zeroTempVEV = [(Higgsfields[i], higgsvevAtZeroTemp[i]) for i in range(len(Higgsfields))]
finiteTempVEV = [
    (Higgsfields[i], higgsVEVAtFiniteTemp[i]) for i in range(len(Higgsfields))
]
fieldsZero = [(x, 0) for x in Higgsfields]

phi1Sq = simplify((Dagger(phi1) * phi1)[0])
phi2Sq = simplify((Dagger(phi2) * phi2)[0])
phi12 = simplify((Dagger(phi1) * phi2)[0])
phi21 = simplify((Dagger(phi2) * phi1)[0])

phi12Sq = simplify(phi12**2)
phi21Sq = simplify(phi21**2)

VHiggsNC = (
    m11sq * phi1Sq
    + m22sq * phi2Sq
    + lambda1 / 2 * phi1Sq**2
    + lambda2 / 2 * phi2Sq**2
    + lambda3 * phi1Sq * phi2Sq
    + lambda4 * phi12 * phi21
)
# We need to expand here otherwise simplify cant handle the simplification well enough
VHiggsHC = expand(
    -m12sq * phi12
    + lambda5 / 2 * phi12**2
    + lambda6 * phi1Sq * phi12
    + lambda7 * phi2Sq * phi12
)
VHiggs = simplify(VHiggsNC + VHiggsHC + conjugate(VHiggsHC))

if im(VHiggs) != 0:
    raise Exception("Higgs potential has an imaginary part with " + str(im(VHiggs)))

# Generate the model
G2HDM = ModelGenerator.ModelGenerator(
    params, dparams, CTTadpoles, Higgsfields, VHiggs, zeroTempVEV, finiteTempVEV
)


# Set Gauge fields
W1, W2, W3, B0 = symbols("W1 W2 W3 B0", real=True)


Dmu = (
    -I * Cg / 2 * (sigma1 * W1 + sigma2 * W2 + sigma3 * W3) - I * Cgs / 2 * sigma0 * B0
)
VGauge = (
    simplify(Dagger(Dmu * phi1) * (Dmu * phi1))[0, 0]
    + simplify(Dagger(Dmu * phi2) * (Dmu * phi2))[0, 0]
)

# Generate Lepton Potentials
NuL = symbols("veL vmuL vtauL", real=True)
ER = symbols("eR muR tauR", real=True)
EL = symbols("eL muL tauL", real=True)
LepBase = NuL + ER + EL


def TypeILeptons():
    ye = sqrt(2) * m_electron / v2
    ymu = sqrt(2) * m_mu / v2
    ytau = sqrt(2) * m_tau / v2

    PiLep = Matrix([[ye, 0, 0], [0, ymu, 0], [0, 0, ytau]])

    VFLep = 0
    for i in range(len(NuL)):
        for j in range(len(ER)):
            VFLep += (NuL[i] * PiLep[i, j] * ER[j]) * phi2[0]

    for i in range(len(EL)):
        for j in range(len(ER)):
            VFLep += (EL[i] * PiLep[i, j] * ER[j]) * phi2[1]

    VFLep = simplify(VFLep)
    return VFLep


def TypeIILeptons():
    ye = sqrt(2) * m_electron / v1
    ymu = sqrt(2) * m_mu / v1
    ytau = sqrt(2) * m_tau / v1

    PiLep = Matrix([[ye, 0, 0], [0, ymu, 0], [0, 0, ytau]])

    VFLep = 0
    for i in range(len(NuL)):
        for j in range(len(ER)):
            VFLep += (NuL[i] * PiLep[i, j] * ER[j]) * phi1[1]

    for i in range(len(EL)):
        for j in range(len(ER)):
            VFLep += (EL[i] * PiLep[i, j] * ER[j]) * phi1[0]

    VFLep = simplify(VFLep)
    return VFLep


if Chosen2HDMType == TwoHDMType.TypeI or Chosen2HDMType == TwoHDMType.Flipped:
    VFLep = TypeILeptons()
else:
    VFLep = TypeIILeptons()


# Generate Quark Potentials
UL = symbols("uL cL tL", real=True)
DL = symbols("dL sL bL", real=True)
UR = symbols("uR cR tR", real=True)
DR = symbols("dR sR bR", real=True)
QuarkBase = UR + DR + UL + DL


def TypeIQuarks():
    yb = sqrt(2) * m_bottom / v2
    yc = sqrt(2) * m_charm / v2
    yd = sqrt(2) * m_down / v2
    ys = sqrt(2) * m_strange / v2
    yt = sqrt(2) * m_top / v2
    yu = sqrt(2) * m_up / v2

    DownCoupling = Matrix([[yd, 0, 0], [0, ys, 0], [0, 0, yb]])
    UpCoupling = Matrix([[yu, 0, 0], [0, yc, 0], [0, 0, yt]])

    ULVector = Matrix([[x] for x in UL])
    DLVector = Matrix([[x] for x in DL])
    URVector = Matrix([[x] for x in UR])
    DRVector = Matrix([[x] for x in DR])

    phiUp = phi2
    phiDown = phi2

    VQuark = ULVector.transpose() * VCKM * DownCoupling * DRVector * phiDown[0]
    VQuark += DLVector.transpose() * DownCoupling * DRVector * phiDown[1]

    VQuark += ULVector.transpose() * UpCoupling * URVector * phiUp[1].conjugate()
    VQuark += (
        -DLVector.transpose()
        * Dagger(VCKM)
        * UpCoupling
        * URVector
        * phiUp[0].conjugate()
    )
    VQuark = simplify(VQuark[0, 0])

    return VQuark


def TypeIIQuarks():
    yb = sqrt(2) * m_bottom / v1
    yd = sqrt(2) * m_down / v1
    ys = sqrt(2) * m_strange / v1

    yt = sqrt(2) * m_top / v2
    yu = sqrt(2) * m_up / v2
    yc = sqrt(2) * m_charm / v2

    DownCoupling = Matrix([[yd, 0, 0], [0, ys, 0], [0, 0, yb]])
    UpCoupling = Matrix([[yu, 0, 0], [0, yc, 0], [0, 0, yt]])

    ULVector = Matrix([[x] for x in UL])
    DLVector = Matrix([[x] for x in DL])
    URVector = Matrix([[x] for x in UR])
    DRVector = Matrix([[x] for x in DR])

    phiUp = phi2
    phiDown = phi1

    VQuark = ULVector.transpose() * VCKM * DownCoupling * DRVector * phiDown[0]
    VQuark += DLVector.transpose() * DownCoupling * DRVector * phiDown[1]

    VQuark += ULVector.transpose() * UpCoupling * URVector * phiUp[1].conjugate()
    VQuark += (
        -DLVector.transpose()
        * Dagger(VCKM)
        * UpCoupling
        * URVector
        * phiUp[0].conjugate()
    )
    VQuark = simplify(VQuark[0, 0])

    return VQuark


if Chosen2HDMType == TwoHDMType.TypeII or Chosen2HDMType == TwoHDMType.Flipped:
    VQuark = TypeIIQuarks()
else:
    VQuark = TypeIQuarks()


# Get the tesnors
G2HDM.setGauge([W1, W2, W3, B0], VGauge)
G2HDM.setLepton(LepBase, VFLep)
G2HDM.setQuark(QuarkBase, VQuark)


def setAdditionalCTEquations():
    # additional equations to define a unique CT solution point
    additionaEquations = []
    additionaEquations.append(dlambda4)
    additionaEquations.append(dRelambda7)
    additionaEquations.append(dRelambda6)
    additionaEquations.append(dImlambda7)
    return additionaEquations


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
        G2HDM.printCTForCPP(setAdditionalCTEquations())
        print("//End CT Calculation")

        print("//Begin CT Order for set_CT_Pot_Par")
        G2HDM.printCTOrder()
        print("//End CT Order")

    if printTensors:
        G2HDM.printModelToCPP()

    if method == "treeSimpl":
        G2HDM.printTreeSimplified()

    if method == "CTSimpl":
        G2HDM.printVCTSimplified()
