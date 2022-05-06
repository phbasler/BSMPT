from sympy import *
from sympy.physics.quantum import *
from sympy.printing.cxxcode import cxxcode

class ModelGenerator:
    _params =None 
    _dparams = None 
    _CTTadpoles = None
    _HiggsFields = None
    _VHiggs =None
    _nHiggs = None
    _VCT = None
    _NablaVCW = None
    _HessianVCW = None
    _replacementLists = {}
    def __init__(self, params, dparams,CTTadpoles,HiggsFields,VHiggs):
        self._params = params
        self._dparams = dparams
        self._CTTadpoles = CTTadpoles
        self._HiggsFields = HiggsFields
        self._VHiggs = VHiggs
        self._nHiggs = len(HiggsFields)
        self._calcCT()
        self._NablaVCW = MatrixSymbol('NCW',self._nHiggs,1)
        self._HessianVCW = MatrixSymbol('HCW',self._nHiggs,self._nHiggs)
        counter=0
        for i in range(self._nHiggs):
            self._replacementLists['NCW[' + str(i) + ']'] = 'NablaWeinberg(' + str(i) + ')'
            for j in range(self._nHiggs):
                self._replacementLists['HCW[' + str(counter) + ']'] = 'HesseWeinberg(' + str(i) + "," + str(j) + ")"
                counter += 1


    def printTensors(self):
        fieldsZero = [(x,0) for x in self._HiggsFields]
        for i in range(self._nHiggs):
            val = diff(self._VHiggs,self._HiggsFields[i])
            val = val.subs(fieldsZero)
            val = simplify(val)
            if val != 0:
                print("Curvature_Higgs_L1[" + str(i) +"] =  " + str(val))

        for i in range(self._nHiggs):
            for j in range(i,self._nHiggs):
                val = diff(self._VHiggs,self._HiggsFields[i],self._HiggsFields[j])
                val = val.subs(fieldsZero)
                val = simplify(val)
                if val != 0:
                    print("Curvature_Higgs_L2[" + str(i) +"]["+ str(j) + "] =  " + str(val))

        for i in range(self._nHiggs):
            for j in range(i,self._nHiggs):
                for k in range(j,self._nHiggs):
                    val = diff(self._VHiggs,self._HiggsFields[i],self._HiggsFields[j],self._HiggsFields[k])
                    val = val.subs(fieldsZero)
                    val = simplify(val)
                    if val != 0:
                        print("Curvature_Higgs_L3[" + str(i) +"]["+ str(j) + "]["+str(k)+"] =  " + str(val))

        for i in range(self._nHiggs):
            for j in range(i,self._nHiggs):
                for k in range(j,self._nHiggs):
                    for l in range(k,self._nHiggs):
                        val = diff(self._VHiggs,self._HiggsFields[i],self._HiggsFields[j],self._HiggsFields[k],self._HiggsFields[l])
                        val = val.subs(fieldsZero)
                        val = simplify(val)
                        if val != 0:
                            print("Curvature_Higgs_L4[" + str(i) +"]["+ str(j) + "]["+str(k)+"]["+str(l)+"] =  " + str(val))

    def _calcCT(self):
        self._VCT=0
        for i in range(len(self._params)):
            self._VCT = self._VCT + self._dparams[i]*diff(self._VHiggs,self._params[i])

        for i in range(len(self._HiggsFields)):
            self._VCT = self._VCT + self._CTTadpoles[i]*self._HiggsFields[i]

        self._VCT = simplify(self._VCT)

    def getVCT(self):
        return self._VCT


    def _CTGetEquations(self,vevRplacement):
        gradient = lambda f, v: Matrix([f]).jacobian(v)
        NablaVCT=simplify(gradient(self._VCT,self._HiggsFields).subs(vevRplacement))
        HessianVCT=simplify(hessian(self._VCT,self._HiggsFields).subs(vevRplacement))
        
        Equations=[]
        for i in range(self._nHiggs):
            val = NablaVCT[i]
            if val != 0:
                Equations.append(val + self._NablaVCW[i])

        for i in range(self._nHiggs):
            for j in range(i,self._nHiggs):
                val = HessianVCT[i,j]
                if val != 0:
                    Equations.append(val+self._HessianVCW[i,j])

        return Equations

    def CTGetMatrixEquation(self,vevRplacement):
        Equations = self._CTGetEquations(vevRplacement)
        dtParamsCombined=self._dparams+[x for x in self._CTTadpoles]
        SysMatrix = zeros(len(Equations),len(dtParamsCombined))
        Target = zeros(len(Equations),1)

        for i in range(len(Equations)):
            Target[i,0] = -Equations[i].subs([(x,0) for x in dtParamsCombined ])
            for j in range(len(dtParamsCombined)):
                SysMatrix[i,j] = diff(Equations[i],dtParamsCombined[j])

        return SysMatrix,Target

    def calcCT(self,vevRplacement):
        Equations = self._CTGetEquations(vevRplacement)
        n,m = self._NablaVCW.shape
        NVCWList = [self._NablaVCW[i,j] for i in range(n) for j in range(m)]
        n,m = self._HessianVCW.shape
        HVCWList = [self._HessianVCW[i,j] for i in range(n) for j in range(i,m)]
        dtParamsCombined=self._dparams+[x for x in self._CTTadpoles]
        extraArgs=NVCWList+HVCWList
        combinedParams=dtParamsCombined+extraArgs
        solution=linsolve(Equations,combinedParams)
        solutionPairs = [(combinedParams[i], solution.args[0][i] ) for i in range(len(dtParamsCombined))]
        identitiesPairs = [(combinedParams[i], solution.args[0][i] ) for i in range(len(dtParamsCombined),len(combinedParams))]
        return solutionPairs, identitiesPairs

    def convertToCPP(self, expr):
        code=cxxcode(expr, standard = 'C++11')
        strToPrint = str(code)
        for key,value in self._replacementLists.items():
            strToPrint=strToPrint.replace(key,value)
        return strToPrint

    def printCTForCPP(self,vevRplacement):
        CTPairs, identities = self.calcCT(vevRplacement)
        for par, val in CTPairs:
            print(self.convertToCPP(par) + " = " + self.convertToCPP(val) + ";")

        