#include <BSMPT/minimizer/Minimizer.h>
#include <map>
#include <vector>
class Compare_N2HDM
{
public:
  using Matrix3D = std::vector<std::vector<std::vector<double>>>;
  using Matrix2D = std::vector<std::vector<double>>;
  Compare_N2HDM();
  Matrix3D CheckTripleCT;
  Matrix3D CheckTripleCW;
  Matrix3D CheckTripleTree;
  std::map<int, BSMPT::Minimizer::EWPTReturnType> EWPTPerSetting;
};
