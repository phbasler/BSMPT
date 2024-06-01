from sympy import (
    symbols,
    Matrix,
    diff,
    simplify,
    Symbol,
    linsolve,
    I,
    hessian,
    zeros,
    expand,
)
from sympy.printing.cxx import cxxcode


class ModelGenerator:
    _params = None
    _dparams = None
    _CTTadpoles = None
    _HiggsFields = None
    _VHiggs = None
    _nHiggs = 0
    _VCT = None
    _NablaVCW = None
    _HessianVCW = None
    _replacementLists = {}
    _VEVAtZeroTemp = None
    _VEVAtFiniteTemp = None
    _GaugeFields = []
    _VGauge = None
    _nGauge = 0
    _LeptonFields = []
    _nLeptons = 0
    _VLep = 0
    _QuarkFields = []
    _nQuarks = 0
    _VQuarks = 0
    _CTParVectorOrder = []

    def __init__(
        self,
        params,
        dparams,
        CTTadpoles,
        HiggsFields,
        VHiggs,
        VEVAtZeroTemp,
        finiteTempVEV,
    ):
        self._params = params
        self._dparams = dparams
        self._CTTadpoles = CTTadpoles
        self._HiggsFields = HiggsFields
        self._VHiggs = VHiggs
        self._nHiggs = len(HiggsFields)
        self._VEVAtZeroTemp = VEVAtZeroTemp
        self._calcVCT()
        self._NablaVCW = Matrix(
            [
                [Symbol("NCW[{},{}]".format(i, j), real=True) for j in range(1)]
                for i in range(self._nHiggs)
            ]
        )
        self._HessianVCW = Matrix(
            [
                [
                    Symbol("HCW[{},{}]".format(i, j), real=True)
                    for j in range(self._nHiggs)
                ]
                for i in range(self._nHiggs)
            ]
        )
        self._TreeLevelTadpoleReplacement = None
        self._calcTreeLevelMinimumConditions()
        self._VEVAtFiniteTemp = finiteTempVEV

        counter = 0
        for i in range(self._nHiggs):
            self._replacementLists["NCW[" + str(i) + "]"] = (
                "NablaWeinberg(" + str(i) + ")"
            )
            for j in range(self._nHiggs):
                self._replacementLists["HCW[" + str(counter) + "]"] = (
                    "HesseWeinberg(" + str(i) + "," + str(j) + ")"
                )
                counter += 1

        for i in range(self._nHiggs):
            self._replacementLists["NCW[" + str(i) + ",0]"] = (
                "NablaWeinberg(" + str(i) + ")"
            )
            for j in range(self._nHiggs):
                self._replacementLists["HCW[" + str(i) + "," + str(j) + "]"] = (
                    "HesseWeinberg(" + str(i) + "," + str(j) + ")"
                )

    def setGauge(self, gaugeFields, VGauge):
        self._GaugeFields = gaugeFields
        self._VGauge = VGauge
        self._nGauge = len(gaugeFields)

    def setLepton(self, leptonFields, VLep):
        self._LeptonFields = leptonFields
        self._nLeptons = len(leptonFields)
        self._VLep = VLep

    def setQuark(self, quarkFields, VQuarks):
        self._nQuarks = len(quarkFields)
        self._VQuarks = VQuarks
        self._QuarkFields = quarkFields

    def printVEVOrder(self):
        vevCounter = 0
        for i in range(self._nHiggs):
            val = self._VEVAtFiniteTemp[i][1]
            if val != 0:
                print("VevOrder.at(" + str(vevCounter) + ") = " + str(i) + ";")
                vevCounter += 1

    def _calcTreeLevelMinimumConditions(self):
        gradient = lambda f, v: Matrix([f]).jacobian(v)
        NablaV = simplify(
            gradient(self._VHiggs, self._HiggsFields).subs(self._VEVAtZeroTemp)
        )
        n, m = NablaV.shape
        Equations = [
            NablaV[i, j] for i in range(n) for j in range(m) if NablaV[i, j] != 0
        ]
        self._TreeLevelTadpoleReplacement = linsolve(Equations, self._params)
        self._TreeLevelTadpoleReplacement = self._TreeLevelTadpoleReplacement.args[0]

    def printTreeLevelMinimumConditions(self):
        for i in range(len(self._params)):
            print(
                self.convertToCPP(self._params[i])
                + " = "
                + self.convertToCPP(self._TreeLevelTadpoleReplacement[i])
                + ";"
            )

    def _printTensorsBasic(self, potential, tensorNamePrefix):
        fieldsZero = [(x, 0) for x in self._HiggsFields]
        for i in range(self._nHiggs):
            val = diff(potential, self._HiggsFields[i])
            val = val.subs(fieldsZero)
            val = simplify(val)
            if val != 0:
                print(
                    tensorNamePrefix
                    + "_L1.at("
                    + str(i)
                    + ") = "
                    + self.convertToCPP(val)
                    + ";"
                )

        for i in range(self._nHiggs):
            for j in range(self._nHiggs):
                val = diff(potential, self._HiggsFields[i], self._HiggsFields[j])
                val = val.subs(fieldsZero)
                val = simplify(val)
                if val != 0:
                    print(
                        tensorNamePrefix
                        + "_L2.at("
                        + str(i)
                        + ").at("
                        + str(j)
                        + ") = "
                        + self.convertToCPP(val)
                        + ";"
                    )

        for i in range(self._nHiggs):
            for j in range(self._nHiggs):
                for k in range(self._nHiggs):
                    val = diff(
                        potential,
                        self._HiggsFields[i],
                        self._HiggsFields[j],
                        self._HiggsFields[k],
                    )
                    val = val.subs(fieldsZero)
                    val = simplify(val)
                    if val != 0:
                        print(
                            tensorNamePrefix
                            + "_L3.at("
                            + str(i)
                            + ").at("
                            + str(j)
                            + ").at("
                            + str(k)
                            + ") = "
                            + self.convertToCPP(val)
                            + ";"
                        )

        for i in range(self._nHiggs):
            for j in range(self._nHiggs):
                for k in range(self._nHiggs):
                    for l in range(self._nHiggs):
                        val = diff(
                            potential,
                            self._HiggsFields[i],
                            self._HiggsFields[j],
                            self._HiggsFields[k],
                            self._HiggsFields[l],
                        )
                        val = val.subs(fieldsZero)
                        val = simplify(val)
                        if val != 0:
                            print(
                                tensorNamePrefix
                                + "_L4.at("
                                + str(i)
                                + ").at("
                                + str(j)
                                + ").at("
                                + str(k)
                                + ").at("
                                + str(l)
                                + ") = "
                                + self.convertToCPP(val)
                                + ";"
                            )

    def printTensors(self):
        self._printTensorsBasic(self._VHiggs, "Curvature_Higgs")

    def printCTTensors(self):
        self._printTensorsBasic(self._VCT, "Curvature_Higgs_CT")

    def _calcVCT(self):
        self._VCT = 0
        for i in range(len(self._params)):
            self._VCT = self._VCT + self._dparams[i] * diff(
                self._VHiggs, self._params[i]
            )

        for i in range(len(self._HiggsFields)):
            newTerm = self._CTTadpoles[i] * self._HiggsFields[i]
            self._VCT = self._VCT + newTerm

        self._VCT = simplify(self._VCT)

    def getVCT(self):
        return self._VCT

    def _CTGetEquations(self):
        gradient = lambda f, v: Matrix([f]).jacobian(v)
        NablaVCT = simplify(
            gradient(self._VCT, self._HiggsFields).subs(self._VEVAtZeroTemp)
        )
        HessianVCT = simplify(
            hessian(self._VCT, self._HiggsFields).subs(self._VEVAtZeroTemp)
        )

        Equations = []
        for i in range(self._nHiggs):
            val = NablaVCT[i]
            if val != 0:
                Equations.append(val + self._NablaVCW[i])

        for i in range(self._nHiggs):
            for j in range(i, self._nHiggs):
                val = HessianVCT[i, j]
                if val != 0:
                    Equations.append(val + self._HessianVCW[i, j])

        return Equations

    def CTGetMatrixEquation(self):
        Equations = self._CTGetEquations()
        dtParamsCombined = self._dparams + [x for x in self._CTTadpoles]
        SysMatrix = zeros(len(Equations), len(dtParamsCombined))
        Target = zeros(len(Equations), 1)

        for i in range(len(Equations)):
            Target[i, 0] = -Equations[i].subs([(x, 0) for x in dtParamsCombined])
            for j in range(len(dtParamsCombined)):
                SysMatrix[i, j] = diff(Equations[i], dtParamsCombined[j])

        return SysMatrix, Target

    def calcCTParams(self, additionalEquations):
        EquationsFromSystem = self._CTGetEquations()
        Equations = EquationsFromSystem + additionalEquations

        n, m = self._NablaVCW.shape
        NVCWList = [self._NablaVCW[i, j] for i in range(n) for j in range(m)]
        n, m = self._HessianVCW.shape
        HVCWList = [self._HessianVCW[i, j] for i in range(n) for j in range(i, m)]
        dtParamsCombined = self._dparams + [x for x in self._CTTadpoles]
        extraArgs = NVCWList + HVCWList
        combinedParams = dtParamsCombined + extraArgs
        solution = linsolve(Equations, combinedParams)
        solutionPairs = [
            (combinedParams[i], solution.args[0][i])
            for i in range(len(dtParamsCombined))
        ]
        identitiesPairs = [
            (combinedParams[i], solution.args[0][i])
            for i in range(len(dtParamsCombined), len(combinedParams))
        ]

        for eq in Equations:
            testedEquation = eq.subs(solutionPairs)
            testedEquation = testedEquation.subs(identitiesPairs)
            testedEquation = simplify(testedEquation)
            if testedEquation != 0:
                raise Exception("No solution for identities found")

        uniqueSolutions = []
        nonUniqueSolutions = []
        for sol in solutionPairs:
            isUnique = True
            for par in self._dparams:
                if sol[1].diff(par) != 0:
                    isUnique = False

            if isUnique:
                uniqueSolutions.append(sol)
            else:
                nonUniqueSolutions.append(sol)

        if len(nonUniqueSolutions) != 0:
            print("Non unique Solutions")
            for lhs, rhs in nonUniqueSolutions:
                print(str(lhs) + " = " + str(rhs))

            raise ValueError(
                "You have non unique solutions in your system. Please define additional equations to choose among them."
            )

        for par, val in solutionPairs:
            self._CTParVectorOrder.append(par)

        return solutionPairs, identitiesPairs

    def convertToCPP(self, expr, assignTo=None):
        II = symbols("II", real=True)
        replExpr = ""
        try:
            replExpr = expr.subs(I, II)
        except AttributeError:
            pass

        custom_functions = {"conjugate": "conj"}

        code = cxxcode(
            replExpr,
            standard="C++17",
            user_functions=custom_functions,
            assign_to=assignTo,
        )
        strToPrint = str(code)
        for key, value in self._replacementLists.items():
            strToPrint = strToPrint.replace(key, value)
        return strToPrint

    def printCTForCPP(self, additionalEquations=[]):
        CTPairs, identities = self.calcCTParams(additionalEquations)

        for par, val in CTPairs:
            print(
                "parCT.push_back("
                + self.convertToCPP(val)
                + "); //"
                + self.convertToCPP(par)
            )

    def printCTOrder(self):
        for i in range(len(self._CTParVectorOrder)):
            par = self._CTParVectorOrder[i]
            print(self.convertToCPP(par) + " = par.at(" + str(i) + ");")

    def printGaugeTensors(self):
        for a in range(self._nGauge):
            for b in range(self._nGauge):
                for i in range(self._nHiggs):
                    for j in range(self._nHiggs):
                        val = diff(
                            self._VGauge,
                            self._GaugeFields[a],
                            self._GaugeFields[b],
                            self._HiggsFields[i],
                            self._HiggsFields[j],
                        )
                        if val != 0:
                            print(
                                "Curvature_Gauge_G2H2.at("
                                + str(a)
                                + ").at("
                                + str(b)
                                + ").at("
                                + str(i)
                                + ").at("
                                + str(j)
                                + ") = "
                                + self.convertToCPP(val)
                                + ";"
                            )

    def printLeptons(self):
        for a in range(self._nLeptons):
            for b in range(self._nLeptons):
                for i in range(self._nHiggs):
                    val = diff(
                        self._VLep,
                        self._LeptonFields[a],
                        self._LeptonFields[b],
                        self._HiggsFields[i],
                    )
                    if val != 0:
                        print(
                            "Curvature_Lepton_F2H1.at("
                            + str(a)
                            + ").at("
                            + str(b)
                            + ").at("
                            + str(i)
                            + ") = "
                            + self.convertToCPP(val)
                            + ";"
                        )

    def printQuarks(self):
        for a in range(self._nQuarks):
            for b in range(self._nQuarks):
                for i in range(self._nHiggs):
                    val = diff(
                        self._VQuarks,
                        self._QuarkFields[a],
                        self._QuarkFields[b],
                        self._HiggsFields[i],
                    )
                    if val != 0:
                        print(
                            "Curvature_Quark_F2H1.at("
                            + str(a)
                            + ").at("
                            + str(b)
                            + ").at("
                            + str(i)
                            + ") = "
                            + self.convertToCPP(val)
                            + ";"
                        )

    def printCtrInfo(self):
        nPar = len(self._params)
        for i in range(len(self._params)):
            lhs = self._params[i]
            rhs = self._TreeLevelTadpoleReplacement[i]
            if lhs != rhs:
                nPar -= 1

        print("nPar = " + str(nPar) + ";")

        lenCT = len(self._dparams) + len(self._CTTadpoles)
        print("nParCT = " + str(lenCT) + ";")
        nVEV = 0
        for i in range(self._nHiggs):
            if self._VEVAtFiniteTemp[i][1] != 0:
                nVEV += 1

        print("nVEV = " + str(nVEV) + ";")

    def printModelToCPP(self):
        print("")
        print("//Begin of constructor info")
        self.printCtrInfo()
        print("//End of constructor info")
        print("")

        print("")
        print("//Begin of VevOrder")
        self.printVEVOrder()
        print("//End of VevOrder")
        print("")

        print("")
        print("//Begin of Tree Level Minimum Relations")
        self.printTreeLevelMinimumConditions()
        print("//End of Tree Level Minimum Relations")
        print("")

        print("")
        print("//Begin of Higgs curvature tensors")
        self.printTensors()
        print("//End of Higgs curvature tensors")
        print("")

        print("")
        print("//Begin of Higgs CT curvature tensors")
        self.printCTTensors()
        print("//End of Higgs CT curvature tensors")
        print("")

        print("")
        print("//Begin of Gauge interaction tensors")
        self.printGaugeTensors()
        print("//End of Gauge interaction tensors")
        print("")

        print("")
        print("//Begin of Lepton interaction tensors")
        self.printLeptons()
        print("//End of Lepton interaction tensors")
        print("")

        print("")
        print("//Begin of Quark interaction tensors")
        self.printQuarks()
        print("//End of Quark interaction tensors")
        print("")

    def printTreeSimplified(self):
        vevCounter = 0
        for i in range(self._nHiggs):
            val = self._VEVAtFiniteTemp[i][1]
            if val != 0:
                print("double " + self.convertToCPP(val) + " = v.at(" + str(i) + ");")
                vevCounter += 1

        VS = simplify(expand(self._VHiggs.subs(self._VEVAtFiniteTemp)))
        print(self.convertToCPP(VS, "res"))

    def printVCTSimplified(self):
        vevCounter = 0
        for i in range(self._nHiggs):
            val = self._VEVAtFiniteTemp[i][1]
            if val != 0:
                print("double " + self.convertToCPP(val) + " = v.at(" + str(i) + ");")
                vevCounter += 1

        VS = simplify(expand(self._VCT.subs(self._VEVAtFiniteTemp)))
        print(self.convertToCPP(VS, "res"))
