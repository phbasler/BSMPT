(* ::Package:: *)

(* ::Section:: *)
(*Disable error message*)


Off[Part::partd]


(* ::Section:: *)
(*Helper functions*)


ToC[x_]:=StringReplace[ToString[CForm[x//HoldForm]],{"Power"->"pow", "vev0" -> "SMConstants.C_vev0", "Sqrt"->"sqrt", "Sin"->"sin", "Cos"->"cos", "Tan"->"tan", "\""->""}]
RemoveDuplicates[list__] := Cases[Tally[list], {x_, 1} :> x]


(* ::Section:: *)
(*Curvatures*)


GenerateCurvature[n_]:= Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[Range[RepHiggsZero//Length],{n}];
CodeLn={};
Table[curvature=(D[VHiggs,Sequence@@higgsbase[[Combination]]]/.RepHiggsZero)//Simplify;
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Higgs_L"<>ToString[n]<>"["<>StringRiffle[Combination-1, "]["]<>"] = "<>ToString[curvature//CForm]<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


(* ::Section:: *)
(*CTs Curvatures*)


GenerateCTsCurvature[n_]:= Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[Range[RepHiggsZero//Length],{n}];
CodeLn={};
Table[curvature=(D[VCT,Sequence@@higgsbase[[Combination]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Higgs_CT_L"<>ToString[n]<>"["<>StringRiffle[Combination-1, "]["]<>"] = "<>ToString[curvature//CForm]<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


(* ::Section:: *)
(*Gauge Curvatures*)


GenerateGaugeCurvature := Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[{Range[Length[GaugeBasis]],Range[Length[GaugeBasis]],Range[Length[higgsbase]],Range[Length[higgsbase]]}];
CodeLn={};
Table[curvature=(D[Vgauge,Sequence@@GaugeBasis[[Combination[[;;2]]]],Sequence@@higgsbase[[Combination[[3;;]]]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Gauge_G2H2["<>StringRiffle[Combination-1, "]["]<>"] = "<>(ToString[(curvature//CForm)]//StringReplace[#,{"Power"->"pow","Subscript(C,g)"->"SMConstants.C_g","Subscript(C,gs)"->"SMConstants.C_gs"}]&)<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


(* ::Section:: *)
(*Leptonic Curvatures*)


GenerateLeptonCurvature := Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[{Range[Length[LepBase]],Range[Length[LepBase]],Range[Length[higgsbase]]}];
CodeLn={};
Table[curvature=(D[VFLep,Sequence@@LepBase[[Combination[[;;2]]]],Sequence@@higgsbase[[Combination[[3;;]]]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Lepton_F2H1["<>StringRiffle[Combination-1, "]["]<>"] = "<>(ToString[(curvature/.RepLepMass//CForm)]//StringReplace[#,{"Power"->"pow","Complex(0,1)"-> "II", "CMass" -> "SMConstants.C_Mass"}]&)<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


(* ::Section:: *)
(*Quark Curvatures*)


GenerateQuarkCurvature := Module[{AllCombinationsN,CodeLn},
AllCombinationsN = Tuples[{Range[Length[baseQuarks]],Range[Length[baseQuarks]],Range[Length[higgsbase]]}];
CodeLn={};
Table[curvature=(D[VF,Sequence@@baseQuarks[[Combination[[;;2]]]],Sequence@@higgsbase[[Combination[[3;;]]]]]/.RepHiggsZero);
If[Not[PossibleZeroQ[curvature]],AppendTo[CodeLn,"Curvature_Quark_F2H1["<>StringRiffle[Combination-1, "]["]<>"] = "<>(ToString[(curvature/.RepMassQuark//CForm)]//StringReplace[#,{"Power"->"pow","Complex(0,1)"-> "II","Complex(0,-1)"-> "-II", "CMass" -> "SMConstants.C_Mass","Conjugate"->"conj"}]&)<>";"]],{Combination,AllCombinationsN}];
Return[CodeLn]];


(* ::Section:: *)
(*Model Implementation*)


ListImplementedModels := Block[{}, models = Import[Nest[ParentDirectory, NotebookDirectory[], 3]<>"/src/models/IncludeAllModels.cpp", "Text"]//StringDelete[#, WhitespaceCharacter]&;
Return[StringSplit[models, {"caseModelIDs::", ":returnstd::make_unique"}][[2;;-2;;2]]]]


InsertCMakeLists[name_] := Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/src/models/CMakeLists.txt";
file = FromCharacterCode[ToCharacterCode[ReadList[filename, "String",NullRecords->True]], "UTF-8"];
If[MemberQ[file,"    ${header_path}/"<> ToString[name] <>".h"] || MemberQ[file,"    ClassPotential"<> ToString[name] <>".cpp"],Print["Model is already in CMakeList.txt"];Return[]];
indexheader = Position[file, "set(header_path \"${BSMPT_SOURCE_DIR}/${suffix}\")"][[1,1]]+1;
indexsrc = Position[file, "set(src"][[1,1]];
indexaddlibrary = Position[file, "add_library(Models STATIC ${header} ${src})"][[1,1]];
before = file[[;;indexheader]];
coreheader = Insert[file[[indexheader+1;; indexsrc-1]],"    ${header_path}/ClassPotential" <> ToString[name] <> ".h",-3];
middlecore = file[[indexsrc;;indexsrc]];
corecpp = Insert[file[[indexsrc+1;;indexaddlibrary-1]],"    ClassPotential" <> ToString[name] <> ".cpp",-3];
after = file[[indexaddlibrary;;]];
file = Join[before, coreheader, middlecore, corecpp, after,{""}];
Export[filename, file, "Table"];]


InsertIncludeAllModelsHeader[name_] := Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/include/BSMPT/utility/ModelIDs.h";
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


InsertIncludeAllModelssrc[name_] := Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/src/models/IncludeAllModels.cpp";
file = FromCharacterCode[ToCharacterCode[ReadList[filename, "String",NullRecords->True]], "UTF-8"];
If[MemberQ[file,"#include <BSMPT/models/ClassPotential" <> ToString[name] <> ".h>"] || MemberQ[file,"  case ModelIDs::" <> ToUpperCase[ToString[name]] <> ":"],Print["Model is already in IncludeAllModels.h"];Return[]];
indexinclude = Position[file, "#include <BSMPT/models/IncludeAllModels.h>"][[1,1]];
indexchoice = Position[file, "  switch (choice)"][[1,1]]+1;
indexinvalidmodel = Position[file, "  default: throw std::runtime_error(\"Invalid model\");"][[1,1]];
include = Insert[file[[;;indexinclude]],"#include <BSMPT/models/ClassPotential" <> ToString[name] <> ".h>",-2];
untilchoice = file[[indexinclude+1;;indexchoice]];
cases = Append[file[[indexchoice+1;;indexinvalidmodel-1]], "  case ModelIDs::" <> ToUpperCase[ToString[name]] <> ":\n    return std::make_unique<Class_Potential_" <> ToString[name] <> ">(smConstants);\n    break;"];
end = file[[indexinvalidmodel;;]];
file = Join[include,untilchoice, cases, end,{""}];
Export[filename, file, "Table"];]


CreateModelHeader[name_,InputParameters_,DepententParameters_,parCT_,higgsvev_] :=Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/include/BSMPT/models/ClassPotential" <> ToString[name] <> ".h";
uppercasename = ToUpperCase[ToString[name]];

file={"#ifndef SRC_CLASSPOTENTIAL"<>uppercasename<>"_H_
#define SRC_CLASSPOTENTIAL"<>uppercasename<>"_H_

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{
class Class_Potential_"<> ToString[name] <>" : public Class_Potential_Origin
{
public:
  Class_Potential_"<>ToString[name]<>"(const ISMConstants &smConstants);
  virtual ~Class_Potential_"<>ToString[name]<>"();
", 
  "  // Initialize input parameters", Sequence@@Table["  double " <> ToString[i] <> " = 0;",{i,InputParameters}],"",
  "  // Initialize dependent parameters", Sequence@@Table["  double " <> ToString[i[[1]]] <> " = 0;",{i,DepententParameters}],"",
  "  // Initialize counter terms", Sequence@@Table["  double " <> ToString[i] <> " = 0;",{i,parCT}],"","
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


CreateModelFile[name_,higgsbase_,higgsvev_,higgsvevFiniteTemp_,VEVList_,par_,InputParameters_,DepententParameters_,CurvatureL1_,CurvatureL2_,CurvatureL3_,CurvatureL4_,GaugeCurvatureL4_,LeptonCurvatureL3_,QuarkCurvatureL3_,parCT_,CTCurvatureL1_,CTCurvatureL2_,CTCurvatureL3_,CTCurvatureL4_,GaugeBasis_,LepBase_,baseQuarks_] :=Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/src/models/ClassPotential" <> ToString[name] <> ".cpp";
uppercasename = ToUpperCase[ToString[name]];

file={"#include \"Eigen/Dense\"
#include \"Eigen/Eigenvalues\"
#include \"Eigen/IterativeLinearSolvers\"
#include <BSMPT/models/SMparam.h> // for SMConstants.C_vev0, SMConstants.C_MassTop, SMConstants.C_g
#include <algorithm> // for max, copy
#include <iomanip>
#include <iostream> // for operator<<, endl, basic_o...
#include <memory>   // for allocator_traits<>::value...
#include <stddef.h> // for std::size_t

#include <BSMPT/models/ClassPotential" <> name <> ".h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{

Class_Potential_" <> name <> "::Class_Potential_" <> name <> "(
    const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model         = ModelID::ModelIDs::" <> uppercasename <> ";

  nPar = " <> ToString[Length[par]] <> ";   // number of parameters in the tree-Level Lagrangian AFTER using
               // tadpole equations
  nParCT = " <> ToString[Length[parCT]] <> "; // number of parameters in the counterterm potential

  nVEV = " <> ToString[Length[Select[higgsvevFiniteTemp,Not[PossibleZeroQ[#]]&]]] <>"; // number of VEVs to minimize the potential

  NHiggs = " <> ToString[Length[higgsbase]] <> "; // number of scalar d.o.f.

  NGauge = " <> ToString[Length[GaugeBasis]] <> "; // number of gauge fields

  NLepton = " <> ToString[Length[LepBase]] <> "; // number of lepton fields

  NQuarks = " <> ToString[Length[baseQuarks]] <> "; // number of quark fields

  VevOrder.resize(nVEV);", Sequence@@Table["  VevOrder[" <> ToString[i-1] <> "] = "<> ToString[VEVList[[i]]] <>"; // " <> ToString[higgsvevFiniteTemp[[VEVList[[i]]+1]]],{i,Length[VEVList]}],
  "
  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_"<>name<>"::~Class_Potential_"<>name<>"()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_"<>name<>"::addLegendCT() const
{
  std::vector<std::string> labels;",
  Sequence@@Table["  labels.push_back(\""<> ToString[i] <>"\");",{i,parCT}],
  "
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_Potential_"<>name<>"::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back(\"T_c\");     // Label for the critical temperature
  labels.push_back(\"v_c\");     // Label for the critical vev
  labels.push_back(\"v_c/T_c\"); // Label for xi_c
  // out += \"VEV order\";",
  Sequence@@Table["  labels.push_back(\"" <> ToString[i] <> "(T_c)\");",{i,Select[higgsvevFiniteTemp,Not[PossibleZeroQ[#]]&]}]
  ,"
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 */
std::vector<std::string>
Class_Potential_"<>name<>"::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  // mass basis, you can identify here your particles",Sequence@@Table["  particles.push_back(\"h_" <> ToString[i] <> "\");",{i,nHiggs}],"
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = i; j < NHiggs; j++)
    {
      for (std::size_t k = j; k < NHiggs; k++)
      {
        labels.push_back(\"Tree_\" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back(\"CT_\" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back(\"CW_\" + particles.at(i) + particles.at(j) +
                         particles.at(k));
      }
    }
  }

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_"<>name<>"::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = \"Your VEV order\"",
  Sequence@@Table["  labels.push_back(\"" <> ToString[i] <> "\");",{i,Select[higgsvevFiniteTemp,Not[PossibleZeroQ[#]]&]}]
  ,"
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_"<>name<>"::ReadAndSet(const std::string &linestr,
                                             std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= "<>ToString[Length[InputParameters]]<>"; k++)
  {
    ss >> tmp;",Sequence@@Table["    if (k == " <> ToString[i] <>")\n      par[" <> ToString[i-1] <> "] = tmp; // " <> ToString[InputParameters[[i]]],{i,Length[InputParameters]}],"
  }

  set_gen(par);
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_"<>name<>"::set_gen(const std::vector<double> &par)
{
", Sequence@@Table["  " <> ToString[InputParameters[[i]]] <>" = par[" <> ToString[i-1] <> "]; ",{i,Length[InputParameters]}],"",Sequence@@Table["  " <> (DepententParameters[[i]][[1]]//ToC) <> " = " <> ToString[DepententParameters[[i]][[2]]//ToC] <> "; ",{i,Length[DepententParameters]}],"
  scale = SMConstants.C_vev0; // renormalisation scale is set to the SM VEV

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // set the vector vevTreeMin. vevTree will then be set by the
  // function MinimizeOrderVEV", Sequence@@Table["  vevTreeMin[" <> ToString[i-1] <> "] = "<> ToString[higgsvev[[VEVList[[i]]+1]]] <>"; // " <> ToString[higgsvevFiniteTemp[[VEVList[[i]]+1]]],{i,Length[VEVList]}],"
  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_Potential_"<>name<>"::set_CT_Pot_Par(const std::vector<double> &par)
{",
  Sequence@@Table["  " <> ToString[parCT[[i]]] <> " = par[" <> ToString[i-1]<>"];",{i,Length[parCT]}],"
  // assign the non-zero entries",
  Sequence@@Table[ "  " <> ToString[i],{i,CTCurvatureL1}] ,"",
  Sequence@@Table[ "  " <> ToString[i],{i,CTCurvatureL2}] ,"",
  Sequence@@Table[ "  " <> ToString[i],{i,CTCurvatureL3}] ,
  Sequence@@Table[ "  " <> ToString[i],{i,CTCurvatureL4}] ,""
  ,"
}

/**
 * console output of all parameters
 */
void Class_Potential_"<>name<>"::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << \"Model = \" << Model << \"\\n\";\n
  ss << \"\\nThe input parameters are : \\n\";",
  Sequence@@Table["  ss << \""<> ToString[i] <> " = \" << " <> ToString[i] <> " << \"\\n\";",{i,InputParameters}],"
  ss << \"\\nThe parameters are : \\n\";",
  Sequence@@Table["  ss << \""<> ToString[i] <> " = \" << " <> ToString[i] <> " << \"\\n\";",{i,par}],"
  ss << \"\\nThe counterterm parameters are : \\n\";",
  Sequence@@Table["  ss << \""<> ToString[i] <> " = \" << " <> ToString[i] <> " << \"\\n\";",{i,parCT}],"
  ss << \"\\nThe scale is given by mu = \" << scale << \" GeV \\n\";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_Potential_"<>name<>"::calc_CT() const
{
  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += \" was called before SetCurvatureArrays()!\\n\";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsdone)
  {
    std::string retmes = __func__;
    retmes += \" was called before CalculatePhysicalCouplings()!\\n\";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  // formulae for the counterterm scheme",Sequence@@Table["  parCT.push_back(" <> ToC[tiCTs[[i]]/.{NCW[x_]->NablaWeinberg[x-1],HCW[x_,y_]->HesseWeinberg[x-1,y-1]}]<> "); //"<>ToString[parCT[[i]]//CForm]<>";",{i,Length[parCT]}],"
  return parCT;
}

// mass basis triple couplings
void Class_Potential_"<>name<>"::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  // new rotation matrix with
  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  std::vector<double> HiggsOrder(NHiggs);

  // example for keeping the mass order
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          for (std::size_t m = 0; m < NHiggs; m++)
          {
            for (std::size_t n = 0; n < NHiggs; n++)
            {
              double RotFac =
                  HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
              TripleHiggsCorrectionsCWPhysical[i][j][k] +=
                  RotFac * GaugeBasis[l][m][n];
              TripleHiggsCorrectionsTreePhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3[l][m][n];
              TripleHiggsCorrectionsCTPhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3_CT[l][m][n];
            }
          }
        }
      }
    }
  }
}

void Class_Potential_"<>name<>"::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // assign the non-zero entries",
  Sequence@@Table[ "  " <> ToString[i],{i,CurvatureL1}] ,
  Sequence@@Table[ "  " <> ToString[i],{i,CurvatureL2}] ,"",
  Sequence@@Table[ "  " <> ToString[i],{i,CurvatureL3}] ,
  Sequence@@Table[ "  " <> ToString[i],{i,CurvatureL4}] ,"",
  Sequence@@Table[ "  " <> ToString[i],{i,GaugeCurvatureL4}] ,
  "
  std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = SMConstants.C_Vud;
  V12 = SMConstants.C_Vus;
  V13 = SMConstants.C_Vub;
  V21 = SMConstants.C_Vcd;
  V22 = SMConstants.C_Vcs;
  V23 = SMConstants.C_Vcb;
  V31 = SMConstants.C_Vtd;
  V32 = SMConstants.C_Vts;
  V33 = SMConstants.C_Vtb;
",
  Sequence@@Table[ "  " <> ToString[i],{i,LeptonCurvatureL3}] ,"",
  Sequence@@Table[ "  " <> ToString[i],{i,QuarkCurvatureL3}] ,"
}

bool Class_Potential_"<>name<>"::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Potential_"<>name<>"::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */
  return false;
}
double
Class_Potential_"<>name<>"::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = " <> ToC[Simplify[VHiggs/.Thread[higgsbase->Table[If[PossibleZeroQ[higgsvevFiniteTemp[[i]]],0, "v[" <> ToString[i-1]<>"]"],{i,Length[higgsvevFiniteTemp]}]]]] <> ";
  
  return res;
}

double Class_Potential_"<>name<>"::VCounterSimplified(
    const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = " <> ToC[Simplify[VCT/.Thread[higgsbase->Table[If[PossibleZeroQ[higgsvevFiniteTemp[[i]]],0, "v[" <> ToString[i-1]<>"]"],{i,Length[higgsvevFiniteTemp]}]]]] <> ";
  
  return res;
}

void Class_Potential_"<>name<>"::Debugging(const std::vector<double> &input,
                                            std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
"};
Export[filename, file, "Table"];]


CreateExamplePoint[name_,InputParameters_,InputParametersExample_] := Block[{},
filename = Nest[ParentDirectory, NotebookDirectory[], 3]<>"/example/" <> name <> "_Input.tsv";
file = {Prepend[InputParameters,Null],{1,Sequence@@InputParametersExample}};
fileTest = FileNames["Test",Nest[ParentDirectory, NotebookDirectory[], 3]<>"/build/*/bin"][[1]]; 
WriteString[$Output,"To run example point \n" <>  fileTest  <> " --model=" <> name <> " --input="<> filename <> " --line=2"];
Export[filename, file, "Table"];]





ImplementModel[
name_,
higgsbase_,
higgsvev_,
higgsvevFiniteTemp_,
VEVList_,par_,
InputParameters_,
InputParametersExample_,
DepententParameters_,
CurvatureL1_,
CurvatureL2_,
CurvatureL3_,
CurvatureL4_,
GaugeCurvatureL4_,
LeptonCurvatureL3_,
QuarkCurvatureL3_,
parCT_,
CTCurvatureL1_,
CTCurvatureL2_,
CTCurvatureL3_,
CTCurvatureL4_,
GaugeBasis_,
LepBase_,
baseQuarks_]:=Block[{},
Print["Inserting model name into CMakeList"];
InsertCMakeLists[name];
Print["Inserting model name into IncludeAllModels.h"];
InsertIncludeAllModelsHeader[name];
Print["Inserting model name into IncludeAllModels.cpp"];
InsertIncludeAllModelssrc[name];
Print["Creating model header"];
CreateModelHeader[name,InputParameters,DepententParameters,parCT,higgsvev];
Print["Creating model implementation .cpp"];
CreateModelFile[name,
higgsbase,
higgsvev,
higgsvevFiniteTemp,
VEVList,
par,
InputParameters,
DepententParameters,
CurvatureL1,
CurvatureL2,
CurvatureL3,
CurvatureL4,
GaugeCurvatureL4,
LeptonCurvatureL3,
QuarkCurvatureL3,
parCT,
CTCurvatureL1,
CTCurvatureL2,
CTCurvatureL3,
CTCurvatureL4,
GaugeBasis,
LepBase,
baseQuarks];
Print["Created example point in BSMPT/example"];
CreateExamplePoint[name,InputParameters,InputParametersExample];
Print[Style["Remember to recompile!",Red]]]


Print[Style["Helper functions successfully imported!",FontSize->20, Directive[Thick]]];
