(* ::Package:: *)

(* ::Chapter:: *)
(*In this section the scalar potential is defined. As an example the N2HDM is shown*)


(*Define higgs fields*)
higgsbase = {\[Rho]1, \[Rho]2, \[Eta]1, \[Eta]2, \[Psi]1, \[Psi]2, \[Zeta]1, \[Zeta]2, \[Zeta]S};


(*Assign vevs at T=0*)
higgsvev = {0, 0, 0, 0, 0, 0, v1, v2, vs};


(*Assign vevs at T != 0*)
higgsvevFiniteTemp = {0, wcb, 0, 0, 0, wcp, w1, w2, ws};


(*Replacement list for the vevs*)
VEVRep = Table[i[[1]]->i[[2]],{i,Transpose[{higgsbase,higgsvev}]}];


(*Replacement list set fields zero*)
RepHiggsZero = Table[i->0,{i,higgsbase}];


(*Define number of Higgses*)
nHiggs = Length[higgsbase];


(*Define parameters of the potential*)
par = {m11Sq, m22Sq, m12Sq, L1, L2, L3, L4, L5, msSq, L6, L7, L8};


(*Define doublets and scalar fields*)
Phi1 = {higgsbase[[1]] + higgsbase[[3]]*I, higgsbase[[7]] + higgsbase[[5]]*I}/Sqrt[2];
Phi2 = {higgsbase[[2]] + higgsbase[[4]]*I, higgsbase[[8]] + higgsbase[[6]]*I}/Sqrt[2];
S = higgsbase[[9]];
(*Print fields*)
Phi1//MatrixForm//TraditionalForm
Phi2//MatrixForm//TraditionalForm
S//MatrixForm//TraditionalForm


(*Write the potential*)
$Assumptions = higgsbase \[Element] Reals;
C11 = Simplify[ConjugateTranspose[Phi1] . Phi1];
C12 = Simplify[ConjugateTranspose[Phi1] . Phi2];
C21 = Simplify[ConjugateTranspose[Phi2] . Phi1];
C22 = Simplify[ConjugateTranspose[Phi2] . Phi2];
CS = Simplify[Conjugate[S]*S];
VHiggs = par[[1]]*C11 + par[[2]]*C22 - par[[3]]*(C12 + C21) + par[[4]]/2*C11^2 + par[[5]]/2*C22^2 + par[[6]]*C11*C22 + par[[7]]*C12*C21 + par[[8]]/2*(C12^2 + C21^2) + par[[9]]/2*CS + par[[10]]/8*CS^2 + par[[11]]/2*C11*CS + par[[12]]/2*C22*CS//Simplify


(*Calculate Tadpoles for Minimisation conditions*)
VHiggsGrad=D[VHiggs,{higgsbase}]/.VEVRep//Simplify;
VHiggsGrad//TableForm


(*Which parameters should be replaced through the minimum conditions?*)
ParToReplace = {m11Sq, m22Sq, msSq};
TadpoleRep = Solve[VHiggsGrad == Table[0,{VHiggsGrad}],ParToReplace];
(DepententParameters = TadpoleRep[[1]])//TableForm


(* ::Chapter:: *)
(*Code*)
(*In this section the model implementation is generated.*)


(* ::Section:: *)
(*MinimizeOrderVEV*)


(*This section calculates the vector VevOrder which is set in the MinimizeOrderVEV function.*)
VEVList={};
Table[If[Not[PossibleZeroQ[higgsvevFiniteTemp[[i]]]],AppendTo[VEVList,i-1]];{i,higgsvevFiniteTemp[[i]]},{i,Length[higgsvevFiniteTemp]}];
Table["VevOrder[" <> ToString[vev-1] <> "] = " <> ToString[VEVList[[vev]]] <> ";" ,{vev,Length[VEVList]}]//TableForm


(* ::Section:: *)
(*Curvatures*)


GenerateCurvature[n_]:= Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[Range[RepHiggsZero//Length],{n}];
CodeLn={};
Table[curvature=(D[VHiggs,Sequence@@higgsbase[[Combination]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Higgs_L"<>ToString[n]<>"["<>StringRiffle[Combination-1, "]["]<>"] = "<>ToString[curvature//CForm]<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


(* ::Section:: *)
(*Calculate Higgs Curvature L1*)


CurvatureL1 = GenerateCurvature[1];
If[CurvatureL1=={},"No entries for Curvature_Higgs_L1",CurvatureL1//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature L2*)


CurvatureL2 = GenerateCurvature[2];
If[CurvatureL2=={},"No entries for Curvature_Higgs_L2",CurvatureL2//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature L3*)


CurvatureL3 = GenerateCurvature[3];
If[CurvatureL3=={},"No entries for Curvature_Higgs_L3",CurvatureL3//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature L4*)


CurvatureL4 = GenerateCurvature[4];
If[CurvatureL4=={},"No entries for Curvature_Higgs_L4",CurvatureL4//TableForm]


(* ::Section:: *)
(*Counterterm potential*)


(*Define counterterms*)
parCT = Join[Table[Symbol["d"<>ToString[p]],{p,par}], Table[Symbol["dT"<>ToString[i]],{i,Length[higgsbase]}]]


(*Define the counterterm potential*)
VCT = ((VHiggs/.Table[par[[i]]->parCT[[i]],{i,Length[par]}]) + parCT[[(par//Length)+1;;]] . higgsbase)


(*Define the sytem of equation to determine the CTs*)
nCount=0;
EqTar={};
Table[nCount+=1; GL[nCount]=D[VCT,higgsbase[[i]]]/.VEVRep; If[PossibleZeroQ[GL[nCount]],nCount-=1,AppendTo[EqTar,-NCW[i]]],{i,Length[higgsbase]}];
Table[nCount+=1; GL[nCount]=D[VCT,higgsbase[[i]],higgsbase[[j]]]/.VEVRep; If[PossibleZeroQ[GL[nCount]],nCount-=1,AppendTo[EqTar,-HCW[i,j]]],{i,Length[higgsbase]},{j,i,Length[higgsbase]}];
EqMatrix = Table[D[GL[i],j],{i,nCount},{j,parCT}];
Print["EqMatrix rank is ", EqMatrix//MatrixRank]
EqMatrix//MatrixForm
Print["EqTar is ", EqTar//TableForm]


(*Define the sytem of equation to determine the CTs*)
SysOriginal = Append[EqMatrix//Transpose,EqTar]//Transpose;
Print["SysOriginal rank is ", SysOriginal//MatrixRank]
SysOriginal//MatrixForm


(*Find the CW relations that help solve the system*)
SysOriginal = Append[EqMatrix//Transpose,EqTar]//Transpose;
Stmp=RowReduce[SysOriginal, ZeroTest -> (! FreeQ[#, HCW[_,_]] || !FreeQ[#,NCW[_]] &)];
CWRelations =Solve[Select[Stmp,#[[;;-2]]==Table[0,{i,21}]&][[All,-1]]==0,-EqTar][[1]];
Print["Relations found between the CW 1st and 2nd derivative\n", CWRelations//TableForm]


(*Reduce the system of equation so that Mathematica can solve it*)
SysOriginal = Append[EqMatrix//Transpose,EqTar]//Transpose;
SysOriginal = (SysOriginal/.CWRelations);
NewMatrix={SysOriginal[[1]]};
nCount=1;
While[MatrixRank[NewMatrix]!=MatrixRank[SysOriginal] && nCount < (Length[SysOriginal] + 1),
nCount+=1;
If[MatrixRank[Append[NewMatrix,SysOriginal[[nCount]]]]>MatrixRank[NewMatrix],AppendTo[NewMatrix,SysOriginal[[nCount]]]]]


(*Determine the CTs*)
NullSpaceofCT = NullSpace[NewMatrix[[All,;;-2]]];
NullSpaceofCT = Sum[NullSpaceofCT[[i]]*Subscript[t, i],{i,Length[NullSpaceofCT]}];
CTs = LinearSolve[NewMatrix[[All,;;-2]],NewMatrix[[All,-1]]];
CTs += NullSpaceofCT;
{parCT,Table["\t\[RightArrow]\t",{i,parCT}],CTs}//Transpose//TableForm


(* ::Subsection:: *)
(*The Subscript[t, i] still need to be fixed *)
(*We chose the following Subscript[t, i] to match the python and Maple implementation, but any other choice is valid .*)


(*Chosen so that dL4 = 0 and dT9 = 0*)
tiChoice = Solve[CTs[[7]]==0 && CTs[[21]]==0,{Subscript[t, 1],Subscript[t, 2]}][[1]]/.CWRelations


(*Fixed CTs*)
tiCTs = CTs/.tiChoice//Simplify;
{parCT,Table["\t\[RightArrow]\t",{i,parCT}],tiCTs}//Transpose//TableForm


(* ::Subsection:: *)
(*Check that the counter terms are correct*)


test = (((((EqMatrix) . parCT)/.(Thread[parCT -> tiCTs]))==(EqTar/.CWRelations))//Simplify);
If[test,Print["Success!"],Print["Something went wrong.\nMake sure the missmatch between the CW derivatives and CT cannot bet solved by a CWRelations substitution rule."];test//TableForm]


(* ::Subsection:: *)
(*Print the CT into C++*)


Table[ToString[parCT[[i]]//CForm] <> " = " <> ToString[CForm[tiCTs[[i]]/.{NCW[x_]->NablaWeinberg[x-1],HCW[x_,y_]->HesseWeinberg[x-1,y-1]}]] <> ";",{i,Length[parCT]}]//TableForm


(* ::Section:: *)
(*CTs Curvatures*)


GenerateCTsCurvature[n_]:= Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[Range[RepHiggsZero//Length],{n}];
CodeLn={};
Table[curvature=(D[VCT,Sequence@@higgsbase[[Combination]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Higgs_CT_L"<>ToString[n]<>"["<>StringRiffle[Combination-1, "]["]<>"] = "<>ToString[curvature//CForm]<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


(* ::Section:: *)
(*Calculate Higgs Curvature CT L1*)


CTCurvatureL1 = GenerateCTsCurvature[1];
If[CTCurvatureL1=={},"No entries for Curvature_Higgs_CT_L1",CTCurvatureL1//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature CT L2*)


CTCurvatureL2 = GenerateCTsCurvature[2];
If[CTCurvatureL2=={},"No entries for Curvature_Higgs_CT_L2",CTCurvatureL2//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature CT L3*)


CTCurvatureL3 = GenerateCTsCurvature[3];
If[CTCurvatureL3=={},"No entries for Curvature_Higgs_CT_L3",CTCurvatureL3//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature CT L4*)


CTCurvatureL4 = GenerateCTsCurvature[4];
If[CTCurvatureL4=={},"No entries for Curvature_Higgs_CT_L4",CTCurvatureL4//TableForm]


(* ::Chapter:: *)
(*Gauge interaction*)
(*This section defines the interaction between gauge bosons and scalars*)
(*First define the basis for the gauge bosons*)


(*Define the gauge fields*)
GaugeBasis={Subscript[W, 1],Subscript[W, 2],Subscript[W, 3],Subscript[B, 0]}


(*Define the covariant derivative*)
Subscript[D, \[Mu]]= -I Subscript[C, g]/2 Sum[Subscript[W, i]PauliMatrix[i],{i,3}] -I Subscript[C, gs]/2 Subscript[B, 0] IdentityMatrix[2];
Subscript[D, \[Mu]]//MatrixForm


(*Define the gauge potential*)
Vgauge = (Conjugate[(Subscript[D, \[Mu]] . Phi1)] . (Subscript[D, \[Mu]] . Phi1) + Conjugate[(Subscript[D, \[Mu]] . Phi2)] . (Subscript[D, \[Mu]] . Phi2))//ComplexExpand//Simplify


(* ::Section:: *)
(*Gauge Curvatures*)


GenerateGaugeCurvature := Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[{Range[Length[GaugeBasis]],Range[Length[GaugeBasis]],Range[Length[higgsbase]],Range[Length[higgsbase]]}];
CodeLn={};
Table[curvature=(D[Vgauge,Sequence@@GaugeBasis[[Combination[[;;2]]]],Sequence@@higgsbase[[Combination[[3;;]]]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Gauge_G2H2["<>StringRiffle[Combination-1, "]["]<>"] = "<>(ToString[(curvature//CForm)]//StringReplace[#,{"Power"->"pow","Subscript(C,g)"->"C_g","Subscript(C,gs)"->"C_gs"}]&)<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


(* ::Section:: *)
(*Calculate Gauge Curvature L4*)


GaugeCurvatureL4 = GenerateGaugeCurvature;
If[GaugeCurvatureL4=={},"No entries for "Curvature_Gauge _G2H2 "",GaugeCurvatureL4//TableForm]


(* ::Chapter:: *)
(*Lepton interaction*)
(*This defines the interactions with Leptons. The N2HDM has several types to avoid FCNC. *)
(*I am showing only Type 2 here as an example, meaning the Leptons and down-type quarks couple to the first doublet while the up-type quark couple to the second doublet.*)


(*First we define the coupling matrix for the leptons*)
PiLep = {{ye, 0, 0}, {0, ymu, 0}, {0, 0, ytau}};
PiLep//MatrixForm


(*Define the components of the lepton sectors*)
NuL = {veL, vmuL, vtauL};
ER = {eR, muR, tauR};
EL = {eL, muL, tauL};
LepBase = {eL, eR, muL, muR, tauL, tauR, veL, vmuL, vtauL};


(*Defining the Lepton potential (=-L_yuk)*)
VFLep = (NuL . PiLep . ER)Phi1[[1]] + (EL . PiLep . ER)Phi1[[2]]


(*Leptonic mass matrix*)
MassLep = D[VFLep,{LepBase,2}]/.VEVRep;
MassLep//MatrixForm


(*Fermionic masses*)
MatrixForm[MLep = Eigenvalues[MassLep]]


(*Calculate the Yukawas*)
RepLepMass = Solve[{MLep[[5]] == CMassElectron, MLep[[7]] == CMassMu, MLep[[9]] == CMassTau}, {ye, ymu, ytau}][[1]]


(* ::Section:: *)
(*Leptonic Curvatures*)


GenerateLeptonCurvature := Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[{Range[Length[LepBase]],Range[Length[LepBase]],Range[Length[higgsbase]]}];
CodeLn={};
Table[curvature=(D[VFLep,Sequence@@LepBase[[Combination[[;;2]]]],Sequence@@higgsbase[[Combination[[3;;]]]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Lepton_F2H1["<>StringRiffle[Combination-1, "]["]<>"] = "<>(ToString[(curvature/.RepLepMass//CForm)]//StringReplace[#,{"Power"->"pow","Complex(0,1)"-> "II", "CMass" -> "C_Mass"}]&)<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


LeptonCurvatureL3 = GenerateLeptonCurvature;
If[LeptonCurvatureL3=={},"No entries for Curvature_Lepton_F2H1",LeptonCurvatureL3//TableForm]


(* ::Chapter:: *)
(*Quark interaction*)
(*As for the Leptons this shows the Type 2 scenario*)


(*First we define the quark doublets*)
UL = {uL, cL, tL};
DL = {dL, sL, bL};
UR = {uR, cR, tR};
DR = {dR, sR, bR};

baseQuarks = {uR, cR, tR, dR, sR, bR, uL, cL, tL, dL, sL, bL};


(*This part is the CKM matrix. If you don't want to include it, set it to unity*)
VCKM = {{V11, V12, V13}, {V21, V22, V23}, {V31, V32, V33}};
VCKM//MatrixForm
CKML = {V11 -> Vud, V12 -> Vus, V13 -> Vub, V21 -> Vcd, V22 -> Vcs, V23 -> Vcb, V31 -> Vtd, V32 -> Vts, V33 -> Vtb};
VCKM/.CKML//MatrixForm


(*Yukawa*)
Dd = {{yd, 0, 0}, {0, ys, 0}, {0, 0, yb}};
DU = {{yu, 0, 0}, {0, yc, 0}, {0, 0, yt}};
Dd//MatrixForm
DU//MatrixForm


VF = (UL . VCKM . Dd . DR)Phi1[[1]] + (DL . Dd . DR)Phi1[[2]] + (UL . DU . UR)Conjugate[Phi2[[2]]] + (DL . ConjugateTranspose[VCKM] . DU . UR)Conjugate[Phi2[[1]]]//Simplify


MQuark = D[VF,{baseQuarks,2}]/.VEVRep;
MQuark//MatrixForm


mQ = Eigenvalues[MQuark];
mQ//MatrixForm


RepMassQuark = Solve[{mQ[[2]] == CMassBottom, mQ[[4]] == CMassCharm, mQ[[6]] == CMassDown, mQ[[8]] ==CMassStrange, mQ[[10]] == CMassTop, mQ[[12]] == CMassUp}, {yb, yc, yd, ys, yt, yu}][[1]];
RepMassQuark//MatrixForm


(* ::Section:: *)
(*Quark Curvatures*)


GenerateQuarkCurvature := Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[{Range[Length[baseQuarks]],Range[Length[baseQuarks]],Range[Length[higgsbase]]}];
CodeLn={};
Table[curvature=(D[VF,Sequence@@baseQuarks[[Combination[[;;2]]]],Sequence@@higgsbase[[Combination[[3;;]]]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Lepton_F2H1["<>StringRiffle[Combination-1, "]["]<>"] = "<>(ToString[(curvature/.RepMassQuark//CForm)]//StringReplace[#,{"Power"->"pow","Complex(0,1)"-> "II","Complex(0,-1)"-> "-II", "CMass" -> "C_Mass","Conjugate"->"conj"}]&)<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


QuarkCurvatureL3 = GenerateQuarkCurvature;
If[QuarkCurvatureL3=={},"No entries for Curvature_Quark_F2H1",QuarkCurvatureL3//TableForm]


ImplementModel[higgsbase_,higgsvev_,higgsvevFiniteTemp_,CurvatureL1_,CurvatureL2_,CurvatureL3_,CurvatureL4_]:=Block[{},1]


ImplementModel[higgsbase, higgsvev, CurvatureL1,CurvatureL2,CurvatureL3,CurvatureL4];


ListImplementedModels := Block[{}, models = Import[Nest[ParentDirectory, NotebookDirectory[], 3]<>"/src/models/IncludeAllModels.cpp", "Text"]//StringDelete[#, WhitespaceCharacter]&;
Return[StringSplit[models, {"caseModelIDs::", ":returnstd::make_unique"}][[2;;-2;;2]]]]
ListImplementedModels


ListImplementedModels


ImplementModel[name_,higgsbase_,higgsvev_,higgsvevFiniteTemp_,CurvatureL1_,CurvatureL2_,CurvatureL3_,CurvatureL4_]:=Block[{},1]


InsertCMakeLists[name_] := Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/src/models/CMakeLists.txt";
file = FromCharacterCode[ToCharacterCode[ReadList[filename, "String",NullRecords->True]], "UTF-8"];
If[MemberQ[file,"    ${header_path}/"<> ToString[name] <>".h"] || MemberQ[file,"    ClassPotential"<> ToString[name] <>".cpp"],Print["Model is already in CMakeList.txt"];Return[]];
indexheader = Position[file, "set(header_path \"${BSMPT_SOURCE_DIR}/include/BSMPT/models\")"][[1,1]]+1;
indexsrc = Position[file, "set(src"][[1,1]];
indexaddlibrary = Position[file, "add_library(Models ${header} ${src})"][[1,1]];
before = file[[;;indexheader]];
coreheader = Insert[file[[indexheader+1;; indexsrc-1]],"    ${header_path}/" <> ToString[name] <> ".h",-3];
middlecore = file[[indexsrc;;indexsrc]];
corecpp = Insert[file[[indexsrc+1;;indexaddlibrary-1]],"    ClassPotential" <> ToString[name] <> ".cpp",-3];
after = file[[indexaddlibrary;;]];
file = Join[before, coreheader, middlecore, corecpp, after,{""}];
Export[filename, file, "Table"];]
InsertCMakeLists["TestModel"]


InsertIncludeAllModelsHeader[name_] := Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/include/BSMPT/models/IncludeAllModels.h";
file = FromCharacterCode[ToCharacterCode[ReadList[filename, "String",NullRecords->True]], "UTF-8"];
If[MemberQ[file,"  " <> ToUpperCase[ToString[name]] <> ","] || MemberQ[file,"    {\"" <> ToLowerCase[ToString[name]] <> "\", ModelIDs::"<> ToUpperCase[ToString[name]] <>"},"],Print["Model is already in IncludeAllModels.h"];Return[]];
indexheader = Position[file, "  // DO NOT EDIT the part below"][[1,1]];
indexunorderedmap = Position[file, "const std::unordered_map<std::string, ModelIDs> ModelNames{"][[1,1]];
indexendunorderedmap = Position[file, "};"][[2,1]];
before = Insert[file[[;;indexheader]],"  " <> ToUpperCase[ToString[name]] <> ",", -3];
middle = file[[indexheader+1;;indexunorderedmap]];
unorderedmap = Append[file[[indexunorderedmap+1;;indexendunorderedmap-1]], "    {\"" <> ToLowerCase[ToString[name]] <> "\", ModelIDs::"<> ToUpperCase[ToString[name]] <>"},"];
end =file[[indexendunorderedmap;;]];
file = Join[before,middle,unorderedmap,end,{""}];
Export[filename, file, "Table"];]
InsertIncludeAllModelsHeader["TestModel"]


InsertIncludeAllModelssrc[name_] := Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/src/models/IncludeAllModels.cpp";
file = FromCharacterCode[ToCharacterCode[ReadList[filename, "String",NullRecords->True]], "UTF-8"];
If[MemberQ[file,"#include <BSMPT/models/ClassPotential" <> ToString[name] <> ".h>"] || MemberQ[file,"  case ModelIDs::" <> ToUpperCase[ToString[name]] <> ":"],Print["Model is already in IncludeAllModels.h"];Return[]];
indexinclude = Position[file, "#include <BSMPT/models/IncludeAllModels.h>"][[1,1]];
indexchoice = Position[file, "  switch (choice)"][[1,1]]+1;
indexinvalidmodel = Position[file, "  default: throw std::runtime_error(\"Invalid model\");"][[1,1]];
include = Insert[file[[;;indexinclude]],"#include <BSMPT/models/ClassPotential" <> ToString[name] <> ".h>",-2];
untilchoice = file[[indexinclude+1;;indexchoice]];
cases = Append[file[[indexchoice+1;;indexinvalidmodel-1]], "  case ModelIDs::" <> ToUpperCase[ToString[name]] <> ":\n    return std::make_unique<Class_" <> ToString[name] <> ">(smConstants);\n    break;"];
end = file[[indexinvalidmodel;;]];
file = Join[include,untilchoice, cases, end,{""}];
Export[filename, file, "Table"];]
InsertIncludeAllModelssrc["TestModel"]


CreateModelHeader[name_,par_,parCT_,higgsvev_] :=Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/include/BSMPT/models/ClassPotential" <> ToString[name] <> ".h";
uppercasename = ToUpperCase[ToString[name]];

file={"#ifndef SRC_CLASSPOTENTIAL"<>uppercasename<>"_H_
#define SRC_CLASSPOTENTIAL"<>uppercasename<>"_H_

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{
class Class_"<> ToString[name] <>" : public Class_Potential_Origin
{
public:
  "<>ToString[name]<>"(const ISMConstants &smConstants);
  virtual ~"<>ToString[name]<>"();
", 
  "  // Initialize variables", Sequence@@Table["  double " <> ToString[i] <> ";",{i,par}],,
  "  // Initialize counter terms", Sequence@@Table["  double " <> ToString[i] <> ";",{i,parCT}],,
  "  // Initialize T = 0 vevs", Sequence@@Table["  double " <> ToString[i] <> ";",{i,Select[higgsvev, Not[PossibleZeroQ[#]]&]}],"
  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
#endif /* SRC_"<>uppercasename<>"_H_ */
"};
Export[filename, file, "Table"];]
CreateModelHeader["TestModel",par,parCT,higgsvev]
