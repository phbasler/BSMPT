// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file minimum tracer class
 */

#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/utility/NumericalDerivatives.h>

using namespace Eigen;

namespace BSMPT
{

std::ostream &operator<<(std::ostream &os, const StatusNLOStability &status)
{
  os << StatusNLOStabilityToString.at(status);
  return os;
}
std::ostream &operator<<(std::ostream &os, const StatusEWSR &status)
{
  os << StatusEWSRToString.at(status);
  return os;
}
std::ostream &operator<<(std::ostream &os, const StatusTracing &status)
{
  os << StatusTracingToString.at(status);
  return os;
}
std::ostream &operator<<(std::ostream &os, const StatusCoexPair &status)
{
  os << StatusCoexPairToString.at(status);
  return os;
}
std::ostream &operator<<(std::ostream &os, const StatusCrit &status)
{
  os << StatusCritToString.at(status);
  return os;
}
std::ostream &operator<<(std::ostream &os, const StatusGW &status)
{
  os << StatusGWToString.at(status);
  return os;
}
std::ostream &operator<<(std::ostream &os, const StatusTemperature &status)
{
  os << StatusTemperatureToString.at(status);
  return os;
}

std::vector<double> MinimumTracer::LocateMinimum(
    const std::vector<double> &guess_In,
    std::function<std::vector<double>(std::vector<double>)> &df,
    std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian,
    const double &error,
    const double &const_multiplier,
    const int &maxiter)
{
  // Save initial guess
  std::vector<double> guess = guess_In;
  // Checks if guess is close enough
  if (L2NormVector(df(guess)) < error) return guess;

  // If not, performs gradient descent until minima with a maximum of
  // "maxiter" iterations
  int dim                       = guess.size(); // Guess dimension
  int i                         = 0;            // Counter
  std::vector<double> new_guess = guess;        // First guess
  std::vector<double> gradient;                 // Gradient
  std::vector<std::vector<double>> Hess;        // HessianNumerical gradient

  for (i = 0; (i < maxiter) && (L2NormVector(df(new_guess)) > error); i++)
  {
    //  Update grad and HessianNumerical :
    gradient = df(new_guess);
    Hess     = Hessian(new_guess);

    // Convert into a Eigen3 // Probably there is a better way
    Eigen::MatrixXd HessMatrix(dim, dim);
    for (int m = 0; m < dim; m++)
    {
      for (int n = 0; n < dim; n++)
      {
        HessMatrix(m, n) = Hess[m][n];
      }
      // Add a small constant to protect us againt numerically unstable negative
      // eigenvalues
      HessMatrix(m, m) += HessianDiagonalShift;
    }

    if (HessMatrix.determinant() != 0) // If HessianNumerical is
    // invertible them do the gradient descent, if not use only the
    // diagonal elements.
    {
      // Convert gradient into vector
      Eigen::VectorXd EigenGradient =
          Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(gradient.data(),
                                                        gradient.size());

      // Calculates the n - th step : x(n + 1) = x(n) + alpha(n)
      Eigen::VectorXd delta =
          HessMatrix.colPivHouseholderQr().solve(-1 * EigenGradient);

      for (int j = 0; j < dim; j++)
      {
        new_guess[j] += delta(j); // Updates guess
      }
    }
    else // If is not invertible use just the gradient
    {
      for (int j = 0; j < dim; j++)
      {
        new_guess[j] -= const_multiplier * gradient[j]; // Updates guess
      }
    }
  }
  return (new_guess);
}

double MinimumTracer::SmallestEigenvalue(
    const std::vector<double> &point,
    const std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian)
{
  std::size_t dim                                  = point.size();
  std::vector<std::vector<double>> current_hessian = Hessian(point);

  Eigen::MatrixXcd mat(dim, dim);
  for (std::size_t i = 0; i < dim; i++)
  {
    mat.col(i) = Eigen::Map<Eigen::VectorXd>(current_hessian[i].data(), dim);
  }

  // Calculate smallest eigenvalue
  double current_min = 1e100;
  for (auto element : mat.eigenvalues())
    current_min = std::min(element.real(), current_min);
  // If the "SmallestEigenvalue" is zero, it can become negative due to
  // numerical errors. To prevent unstable behaviour we add a small constant.
  return current_min + 1e-7;
}

std::vector<double>
MinimumTracer::FindZeroSmallestEigenvalue(std::vector<double> point_1,
                                          double T_1,
                                          std::vector<double> point_2,
                                          double T_2)
{
  double eps = 0.1;
  double ev_1, ev_2, ev_m, T_m; // Eigenvalues of phases and middle temperature
  int dim = this->modelPointer->get_nVEV();
  std::vector<double> point_m;
  std::function<double(std::vector<double>)> V_1, V_2, V_m;
  std::function<std::vector<double>(std::vector<double>)> dV_1, dV_2, dV_m;
  std::function<std::vector<std::vector<double>>(std::vector<double>)>
      Hessian_1, Hessian_2, Hessian_m;

  // Define potential 1
  V_1 = [&](std::vector<double> vev)
  {
    // Potential wrapper
    std::vector<double> res = this->modelPointer->MinimizeOrderVEV(vev);
    return this->modelPointer->VEff(res, T_1) / (1 + T_1 * T_1);
  };
  dV_1      = [&](auto const &arg) { return NablaNumerical(arg, V_1, eps); };
  Hessian_1 = [&](auto const &arg) { return HessianNumerical(arg, V_1, eps); };

  // Define potential 2
  V_2 = [&](std::vector<double> vev)
  {
    // Potential wrapper
    std::vector<double> res = this->modelPointer->MinimizeOrderVEV(vev);
    return this->modelPointer->VEff(res, T_2) / (1 + T_2 * T_2);
  };
  dV_2      = [=](auto const &arg) { return NablaNumerical(arg, V_2, eps); };
  Hessian_2 = [=](auto const &arg) { return HessianNumerical(arg, V_2, eps); };

  // Initial guess for middle point
  point_m = point_1;

  // Calculate smallest EV
  ev_1 = SmallestEigenvalue(point_1, Hessian_1);
  ev_2 = SmallestEigenvalue(point_2, Hessian_2);

  if (ev_1 * ev_2 > 0)
  {
    // Safety check. We should never be here.
    std::stringstream ss;
    ss << "They have the same sign, no zero to be found!\t" << ev_1 << "\t"
       << ev_2 << "\n";
    Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());
    point_1.push_back(-1);
    return point_1;
  }

  // This is to provent instabilities and make sure we do not leave the
  // phase
  double AllowedMaxDistance = fmax(1, L2NormVector(point_1 - point_2) / dim);
  std::vector<double> CandidatePoint = point_1;
  double CandidateTemperature        = T_1;

  while (abs(T_1 / T_2 - 1) > 1e-8)
  {
    T_m = (T_1 + T_2) / 2.;
    // Define potential in the middle
    V_m = [&](std::vector<double> vev)
    {
      // Potential wrapper
      std::vector<double> res = this->modelPointer->MinimizeOrderVEV(vev);
      return this->modelPointer->VEff(res, T_m) / (1 + T_m * T_m);
    };
    dV_m      = [=](auto const &arg) { return NablaNumerical(arg, V_m, eps); };
    Hessian_m = [=](auto const &arg)
    { return HessianNumerical(arg, V_m, eps); };
    point_m =
        LocateMinimum(point_m, dV_m, Hessian_m, 1e-3 * dim / (1 + T_m * T_m));
    ev_m = SmallestEigenvalue(point_m, Hessian_m);

    // Check if there is numerical instabilities so reduce the searching area.
    if (L2NormVector(point_m - point_2) / dim > AllowedMaxDistance or
        L2NormVector(point_1 - point_m) / dim > AllowedMaxDistance)
    {
      T_2 = T_m;
    }
    else if (ev_1 * ev_m > 0)
    {
      // This point has the correct eigenvalue and is close.
      CandidatePoint       = point_m;
      CandidateTemperature = T_m;
    }

    if (ev_1 * ev_m > 0)
    {
      T_1 = T_m;
    }
    else if (ev_2 * ev_m > 0)
    {
      T_2 = T_m;
    }
    else
    {
      // Here ev_m == 0
      break;
    }
  }
  CandidatePoint.push_back(CandidateTemperature);
  return CandidatePoint;
}

std::vector<Minimum>
MinimumTracer::TrackPhase(double &globMinEndT,
                          const std::vector<double> &point_In,
                          const double &currentT_In,
                          const double &finalT,
                          const double &dT_In,
                          const bool &output,
                          const bool &unprotected)
{
  // Test phase tracker
  int dim         = this->modelPointer->get_nVEV();
  int IsInMin     = 0;
  double initialT = currentT_In;
  double currentT = currentT_In;
  double dT       = dT_In;
  double initialdT;
  double eps = 0.1;
  double LengthGradient, PotentialDifference, Distance;
  std::function<std::vector<double>(std::vector<double>)> dV;
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian;
  std::vector<double> point     = point_In;
  std::vector<double> new_point = point_In;
  std::vector<double> zeroTemp;
  std::stringstream ss;
  std::vector<Minimum> MinimumList;
  Minimum newMinimum;

  bool old_min_is_global = true;

  if ((finalT - currentT) / dT < 0)
  {
    // Step has the wrong sign
    dT *= -1;
  }
  initialdT = dT;

  // Reduce the VEV into the same sector
  ReduceVEV(point);

  // Remove flat directions
  ConvertToNonFlatDirections(point);

  ss << "Tracking phase from T = " << currentT << " to T = " << finalT
     << " | Starting minimum at = " << point << "\n";
  while ((finalT - currentT) / dT >= 0)
  {
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      std::vector<double> res = this->modelPointer->MinimizeOrderVEV(vev);
      return this->modelPointer->VEff(res, currentT) /
             (1 + currentT * currentT);
    };
    dV      = [=](auto const &arg) { return NablaNumerical(arg, V, eps); };
    Hessian = [=](auto const &arg) { return HessianNumerical(arg, V, eps); };

    // Locate the minimum
    new_point =
        LocateMinimum(point, dV, Hessian, 1e-4 * GradientThreshold * dim);

    // Reduce the VEV into the same sector
    ReduceVEV(new_point);

    // Remove flat directions
    ConvertToNonFlatDirections(new_point);

    // Calculate the length of the gradient in the normal potential divided
    // by the dimension of the VEV space
    LengthGradient =
        L2NormVector(dV(new_point)) / dim; // (1 + currentT * currentT) *
    // Compare minimum and previous iteration
    Distance = L2NormVector(new_point - point);
    // Compute difference in energy between both minimum
    PotentialDifference = V(new_point) - V(point);

    // If the minimum it n GeV away then we consider it a new phase
    double ThresholdDistance = (double)dim;
    if (abs(dT) < 1e-5)
    {
      break;
    }
    else if (LengthGradient > GradientThreshold or
             Distance > ThresholdDistance or
             std::any_of(new_point.begin(),
                         new_point.end(),
                         [](double v) { return isnan(v); }))
    {
      // Algorithm flew away. Going back one iteration.
      // Decrease temperature step size
      // Not a stationary point anymore. Going back one iteration.
      // Decrease temperature step size

      ss << "\n\033[1;95m.-> | T = " << currentT
         << " |Grad|/dim =  " << LengthGradient << " | Distance = " << Distance
         << " | |dphi/dT| = " << abs(Distance / dT)
         << " | deltaV = " << PotentialDifference
         << " | SEV = " << SmallestEigenvalue(new_point, Hessian)
         << " | <-\033[0m\n";

      if (initialT == currentT)
      {
        ss << "Could not locate the starting minimum at T = " << initialT
           << " GeV.";
        ss << " |Grad|/dim =  " << L2NormVector(dV(point)) / dim
           << " | Distance = " << Distance
           << " | |dphi/dT| = " << abs(Distance / dT)
           << " | deltaV = " << PotentialDifference
           << " | SEV = " << SmallestEigenvalue(point, Hessian);
        if (output) Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());
        return MinimumList;
      }

      currentT -= dT;
      dT /= 10.;
    }
    else if (unprotected)
    {
      newMinimum.point     = new_point;
      newMinimum.temp      = currentT;
      newMinimum.potential = V(new_point) * (1 + currentT * currentT);
      MinimumList.push_back(newMinimum);
      point = new_point;
      // Sucess minimum!
      ss << "\033[1;32m.\033[0m";
      dT *=
          .5 * ThresholdDistance / Distance; // Try to predict the best stepsize
      if (abs(initialdT) <= abs(dT)) dT = initialdT;
    }
    else
    {
      // It is a nearby stationary point!
      if (SmallestEigenvalue(new_point, Hessian) < 0)
      {
        if (IsInMin == -1)
        {
          zeroTemp = FindZeroSmallestEigenvalue(
              point, currentT - dT, new_point, currentT);
          if (zeroTemp.back() > 0)
          {
            currentT = zeroTemp.back();
            zeroTemp.pop_back();
            // Reduce the VEV into the same sector
            ReduceVEV(zeroTemp);
            // Remove flat directions
            ConvertToNonFlatDirections(zeroTemp);
            newMinimum.point     = zeroTemp;
            newMinimum.temp      = currentT;
            newMinimum.potential = V(zeroTemp) * (1 + currentT * currentT);

            if (old_min_is_global) // check if newMinimum is
                                   // still global minimum
            {
              IsGlobMin(newMinimum);
              if (!newMinimum.is_glob_min)
              {
                if (output)
                  Logger::Write(
                      LoggingLevel::MinTracerDetailed,
                      "Phase no longer coincides with global minimum at " +
                          std::to_string(newMinimum.temp));
                globMinEndT       = newMinimum.temp;
                old_min_is_global = false;
              }
            }

            MinimumList.push_back(newMinimum);
          }
          // End tracking
          break;
        }
        else
        {
          ss << "Calculation of phase tracker failed. T  = " << currentT
             << " GeV\t|\t Final T = " << finalT << " GeV\t|\t"
             << LengthGradient << "\t|\t"
             << SmallestEigenvalue(new_point, Hessian) << "\t";
          // Sucess saddle point!
          ss << "\033[1;31m.\033[0m";
          if (output) Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());
          ss.str(std::string());
          return MinimumList; // return if starts in a saddle point
        }
        IsInMin = 1;
      }
      else
      {
        if (IsInMin == 1)
        {
          // We should never be here. We only track minimum
          break;
        }
        else
        {
          newMinimum.point     = new_point;
          newMinimum.temp      = currentT;
          newMinimum.potential = V(new_point) * (1 + currentT * currentT);

          if (old_min_is_global) // check if newMinimum is
                                 // still global minimum
          {
            IsGlobMin(newMinimum);
            if (!newMinimum.is_glob_min)
            {
              if (output)
                Logger::Write(
                    LoggingLevel::MinTracerDetailed,
                    "Phase no longer coincides with global minimum at " +
                        std::to_string(newMinimum.temp));
              globMinEndT       = newMinimum.temp;
              old_min_is_global = false;
            }
          }

          MinimumList.push_back(newMinimum);
          // Sucess minimum!
          ss << "\033[1;32m.\033[0m";
          dT *= .5 * ThresholdDistance /
                Distance; // Try to predict the best stepsize
          if (abs(initialdT) <= abs(dT)) dT = initialdT;
        }
        IsInMin = -1;
      }
      point = new_point;
    }

    // Make sure that or step is not bigger than it should be and we overshot
    // the next minimum
    if (currentT == finalT) break;
    bool SafeStep = abs(finalT - currentT) > abs(dT);
    if (SafeStep)
      currentT += dT;
    else
      currentT = finalT;
  }
  if (output) Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());
  if (output)
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "T = " + std::to_string(currentT) + " GeV");
  if (old_min_is_global) globMinEndT = currentT;
  return MinimumList;
}

std::vector<Minimum>
MinimumTracer::TrackPhase(const std::vector<double> &point_In,
                          const double &currentT_In,
                          const double &finalT,
                          const double &dT_In,
                          const bool &output,
                          const bool &unprotected)
{
  // Test phase tracker
  int dim         = this->modelPointer->get_nVEV();
  int IsInMin     = 0;
  double initialT = currentT_In;
  double currentT = currentT_In;
  double dT       = dT_In;
  double initialdT;
  double eps = 0.1;
  double LengthGradient, PotentialDifference, Distance;
  std::function<std::vector<double>(std::vector<double>)> dV;
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian;
  std::vector<double> point     = point_In;
  std::vector<double> new_point = point;
  std::vector<double> zeroTemp;
  std::stringstream ss;
  std::vector<Minimum> MinimumList;
  Minimum newMinimum;

  if ((finalT - currentT) / dT < 0)
  {
    // Step has the wrong sign
    dT *= -1;
  }
  initialdT = dT;

  // Reduce the VEV into the same sector
  ReduceVEV(point);

  // Remove flat directions
  ConvertToNonFlatDirections(point);

  ss << "Tracking phase from T = " << currentT << " to T = " << finalT
     << " | Starting minimum at = " << point << "\n";
  while ((finalT - currentT) / dT >= 0)
  {
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      std::vector<double> res = this->modelPointer->MinimizeOrderVEV(vev);
      return this->modelPointer->VEff(res, currentT) /
             (1 + currentT * currentT);
    };
    dV      = [=](auto const &arg) { return NablaNumerical(arg, V, eps); };
    Hessian = [=](auto const &arg) { return HessianNumerical(arg, V, eps); };

    // Locate the minimum
    new_point =
        LocateMinimum(point, dV, Hessian, 1e-4 * GradientThreshold * dim);

    // Reduce the VEV into the same sector
    ReduceVEV(new_point);

    // Remove flat directions
    ConvertToNonFlatDirections(new_point);

    // Calculate the length of the gradient in the normal potential divided
    // by the dimension of the VEV space
    LengthGradient =
        L2NormVector(dV(new_point)) / dim; // (1 + currentT * currentT) *
    // Compare minimum and previous iteration
    Distance = L2NormVector(new_point - point);
    // Compute difference in energy between both minimum
    PotentialDifference = V(new_point) - V(point);

    // If the minimum it n GeV away then we consider it a new phase
    double ThresholdDistance = (double)dim;
    if (abs(dT) < 1e-5)
    {
      break;
    }
    else if (LengthGradient > GradientThreshold or
             Distance > ThresholdDistance or
             std::any_of(new_point.begin(),
                         new_point.end(),
                         [](double v) { return isnan(v); }))
    {
      // Algorithm flew away. Going back one iteration.
      // Decrease temperature step size
      // Not a stationary point anymore. Going back one iteration.
      // Decrease temperature step size

      ss << "\n\033[1;95m.-> | T = " << currentT
         << " |Grad|/dim =  " << LengthGradient << " | Distance = " << Distance
         << " | |dphi/dT| = " << abs(Distance / dT)
         << " | deltaV = " << PotentialDifference
         << " | SEV = " << SmallestEigenvalue(new_point, Hessian)
         << " | <-\033[0m\n";

      if (initialT == currentT)
      {
        ss << "Could not locate the starting minimum at T = " << initialT
           << " GeV.";
        ss << " |Grad|/dim =  " << L2NormVector(dV(point)) / dim
           << " | Distance = " << Distance
           << " | |dphi/dT| = " << abs(Distance / dT)
           << " | deltaV = " << PotentialDifference
           << " | SEV = " << SmallestEigenvalue(point, Hessian);
        if (output) Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());
        return MinimumList;
      }

      currentT -= dT;
      dT /= 10.;
    }
    else if (unprotected)
    {
      newMinimum.point     = new_point;
      newMinimum.temp      = currentT;
      newMinimum.potential = V(new_point) * (1 + currentT * currentT);
      MinimumList.push_back(newMinimum);
      point = new_point;
      // Sucess minimum!
      ss << "\033[1;32m.\033[0m";
      dT *=
          .5 * ThresholdDistance / Distance; // Try to predict the best stepsize
      if (abs(initialdT) <= abs(dT)) dT = initialdT;
    }
    else
    {
      // It is a nearby stationary point!
      if (SmallestEigenvalue(new_point, Hessian) < 0)
      {
        if (IsInMin == -1)
        {
          zeroTemp = FindZeroSmallestEigenvalue(
              point, currentT - dT, new_point, currentT);
          if (zeroTemp.back() > 0)
          {
            currentT = zeroTemp.back();
            zeroTemp.pop_back();
            // Reduce the VEV into the same sector
            ReduceVEV(zeroTemp);
            // Remove flat directions
            ConvertToNonFlatDirections(zeroTemp);
            newMinimum.point     = zeroTemp;
            newMinimum.temp      = currentT;
            newMinimum.potential = V(zeroTemp) * (1 + currentT * currentT);
            MinimumList.push_back(newMinimum);
          }
          // End tracking
          break;
        }
        else
        {
          ss << "Calculation of Tc or phase tracker failed. T  = " << currentT
             << " GeV\t|\t Final T = " << finalT << " GeV\t|\t"
             << LengthGradient << "\t|\t"
             << SmallestEigenvalue(new_point, Hessian) << "\t";
          // Sucess saddle point!
          ss << "\033[1;31m.\033[0m";
          if (output) Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());
          ss.str(std::string());
          return MinimumList; // return if starts in a saddle point
        }
        IsInMin = 1;
      }
      else
      {
        if (IsInMin == 1)
        {
          // We should never be here. We only track minimum
          break;
        }
        else
        {
          newMinimum.point     = new_point;
          newMinimum.temp      = currentT;
          newMinimum.potential = V(new_point) * (1 + currentT * currentT);
          MinimumList.push_back(newMinimum);
          // Sucess minimum!
          ss << "\033[1;32m.\033[0m";
          dT *= .5 * ThresholdDistance /
                Distance; // Try to predict the best stepsize
          if (abs(initialdT) <= abs(dT)) dT = initialdT;
        }
        IsInMin = -1;
      }
      point = new_point;
    }

    // Make sure that or step is not bigger than it should be and we overshot
    // the next minimum
    if (currentT == finalT) break;
    bool SafeStep = abs(finalT - currentT) > abs(dT);
    if (SafeStep)
      currentT += dT;
    else
      currentT = finalT;
  }
  if (output) Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());
  if (output)
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "T = " + std::to_string(currentT) + " GeV");
  return MinimumList;
}

void MinimumTracer::ReduceVEV(std::vector<double> &vev)
{
  // Saveguard if GroupElements is not populated
  if (GroupElements.size() == 0) return;
  int MaximumMeasure = -1;
  char *ptr;
  std::string BinaryNumber;
  for (const auto &GroupElement : GroupElements)
  {
    // Clean buffer
    BinaryNumber.clear();
    // Feed buffer
    for (size_t i = 0; i < vev.size(); i++)
      BinaryNumber.append(std::to_string(vev[i] * GroupElement(i, i) >=
                                         1e-5)); // Heaviside function

    // Converto to decimal and compare
    if (std::strtol(BinaryNumber.c_str(), &ptr, 2) > MaximumMeasure)
    {
      MaximumMeasure = std::strtol(BinaryNumber.c_str(), &ptr, 2);
      for (size_t i = 0; i < vev.size(); i++)
        vev[i] *= GroupElement(i, i); // Transform vev
    }

    // Reached maximum measure
    if (MaximumMeasure == pow(2, vev.size()) - 1) return;
  }
}

void MinimumTracer::ReduceVEV(Minimum &min)
{
  ReduceVEV(min.point);
}

const std::vector<std::vector<double>>
MinimumTracer::WarpPath(const std::vector<std::vector<double>> &path,
                        const std::vector<double> &T1,
                        const std::vector<double> &F1,
                        const std::vector<double> &T2,
                        const std::vector<double> &F2)
{
  std::vector<std::vector<double>> r;

  for (std::size_t i = 0; i < path.size(); i++)
  {
    std::vector<double> temp = path[i];
    for (std::size_t d = 0; d < temp.size(); d++)
    {
      // warp the path
      if ((F1[d] - T1[d]) == 0)
        temp[d] = T2[d] + (temp[d] - T1[d]);
      else
        temp[d] = T2[d] + (F2[d] - T2[d]) * (temp[d] - T1[d]) / (F1[d] - T1[d]);
    }
    r.push_back(temp);
  }
  return r;
}

MinimumTracer::MinimumTracer()
{
}

MinimumTracer::MinimumTracer(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const int &WhichMinimizer_in,
    const bool &UseMultithreading_in)
{
  modelPointer      = pointer_in;
  WhichMinimizer    = WhichMinimizer_in;
  UseMultithreading = UseMultithreading_in;

  FindDiscreteSymmetries();
  FindFlatDirections();
}

void MinimumTracer::FindFlatDirections()
{
  // The number 2, 100, 200 were choosen arbitrarily as an example of a S0(3)
  // rotation
  auto nvev      = this->modelPointer->get_nVEV();
  auto vev_order = this->modelPointer->Get_VevOrder();

  std::vector<double> point = std::vector<double>(nvev, 1);

  // initializsation of vector NonFlatDirections
  NonFlatDirections = std::vector<int>(nvev, 1);

  double res_1, res_2;

  // 3-dimensional flat directions (SO(3))
  for (std::size_t i = 0; i < nvev; i++)
  {
    for (std::size_t j = i + 1; j < nvev; j++)
    {
      for (std::size_t k = j; k < nvev; k++)
      {
        point.at(i) = 2;
        point.at(j) = 100;
        point.at(k) = 200;
        res_1       = this->modelPointer->VEff(
            this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
        point.at(i) = 2;
        point.at(j) = 200;
        point.at(k) = 100;
        res_2       = this->modelPointer->VEff(
            this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
        if (almost_the_same(res_1, res_2, 1e-8))
        {
          point.at(i) = 100;
          point.at(j) = 2;
          point.at(k) = 200;
          res_1       = this->modelPointer->VEff(
              this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
          if (almost_the_same(res_1, res_2, 1e-8))
          {
            point.at(i) = 100;
            point.at(j) = 200;
            point.at(k) = 2;
            res_2       = this->modelPointer->VEff(
                this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
            if (almost_the_same(res_1, res_2, 1e-8))
            {
              point.at(i) = 200;
              point.at(j) = 2;
              point.at(k) = 100;
              res_1       = this->modelPointer->VEff(
                  this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
              if (almost_the_same(res_1, res_2, 1e-8))
              {
                point.at(i) = 200;
                point.at(j) = 100;
                point.at(k) = 2;
                res_2       = this->modelPointer->VEff(
                    this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
                if (almost_the_same(res_1, res_2, 1e-8))
                {
                  NonFlatDirections.at(j) = 0; // choose one point in sphere
                  NonFlatDirections.at(k) = 0; // choose one point in sphere
                  flat_dirs_found         = true;
                  flat_3D_dirs.push_back(std::vector<std::size_t>({i, j, k}));
                }
              }
            }
          }
        }

        point.at(i) = 1;
        point.at(j) = 1;
        point.at(k) = 1;
      }
    }
  }

  // 2-dimensional flat directions (SO(2))
  // The number 2, 100 were choosen arbitrarily as an example of a S0(2)
  // rotation
  for (std::size_t i = 0; i < nvev; i++)
  {
    for (std::size_t j = i + 1; j < nvev; j++)
    {
      point.at(i) = 2;
      point.at(j) = 100;
      res_1       = this->modelPointer->VEff(
          this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
      point.at(i) = 100;
      point.at(j) = 2;
      res_2       = this->modelPointer->VEff(
          this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
      point.at(i) = 1;
      point.at(j) = 1;

      if (almost_the_same(res_1, res_2, 1e-8))
      {
        NonFlatDirections.at(j) = 0; // choose one point in circle
        flat_dirs_found         = true;
        flat_2D_dirs.push_back(std::vector<std::size_t>({i, j}));
      }
    }
  }

  // 1-dimensional flat directions
  for (std::size_t i = 0; i < nvev; i++)
  {
    point.at(i) = 2;
    res_1       = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
    point.at(i) = 100;
    res_2       = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(point), 0, 0, 0);
    point.at(i) = 1;

    if (almost_the_same(res_1, res_2, 1e-8))
    {
      NonFlatDirections.at(i) = 0; // choose one field point in well
      flat_dirs_found         = true;
      flat_1D_dirs.push_back(i);
    }
  }

  // todo: higher-dimensional flat directions

  if (flat_dirs_found)
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Flat directions encountered!");
    if (flat_1D_dirs.size() > 0)
    {
      for (auto ind : flat_1D_dirs)
      {
        Logger::Write(LoggingLevel::MinTracerDetailed,
                      "Field index " + std::to_string(ind) +
                          " identified as 1D flat direction.");
      }
    }
    else if (flat_2D_dirs.size() > 0)
    {
      for (auto inds : flat_2D_dirs)
      {
        Logger::Write(LoggingLevel::MinTracerDetailed,
                      "Field index (" + vec_to_string(inds) +
                          ") identified as 2D flat direction.");
      }
    }
    else if (flat_3D_dirs.size() > 0)
    {
      for (auto inds : flat_3D_dirs)
      {
        Logger::Write(LoggingLevel::MinTracerDetailed,
                      "Field index (" + vec_to_string(inds) +
                          ") identified as 3D flat direction.");
      }
    }
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Choosing minimal set of field "
                  "directions as: " +
                      vec_to_string(NonFlatDirections));
  }
  // Remove duplicates
  flat_1D_dirs.erase(unique(flat_1D_dirs.begin(), flat_1D_dirs.end()),
                     flat_1D_dirs.end());
  flat_2D_dirs.erase(unique(flat_2D_dirs.begin(), flat_2D_dirs.end()),
                     flat_2D_dirs.end());
  flat_3D_dirs.erase(unique(flat_3D_dirs.begin(), flat_3D_dirs.end()),
                     flat_3D_dirs.end());
  return;
}

void MinimumTracer::ConvertToNonFlatDirections(std::vector<double> &point)
{
  if (flat_dirs_found) // flat directions in point
  {
    for (const auto &pair : flat_2D_dirs)
    {
      auto a               = point.at(pair.at(0));
      auto b               = point.at(pair.at(1));
      auto c               = std::sqrt(a * a + b * b);
      point.at(pair.at(0)) = c;
      point.at(pair.at(1)) = 0;
    }
  }

  return;
}

void MinimumTracer::FindDiscreteSymmetries()
{
  const double GroupElementslMaximumRelativeError =
      1e-8; // Maximum value of |V/GroupElements(V)-1|
  const size_t NumberOfGroupElementsChecks =
      10; // Number of times GroupElements is checked
  const size_t dim = modelPointer->get_nVEV(); // Number of VEVs

  std::function<double(Eigen::VectorXd)> V = [&](Eigen::VectorXd vev)
  {
    // Potential wrapper at T=0 for tree-level potential
    std::vector<double> res = this->modelPointer->MinimizeOrderVEV(
        std::vector<double>(vev.data(), vev.data() + vev.size()));
    return this->modelPointer->VEff(res, 0, 0, 0);
  };

  // Generate random VEV
  srand(time(NULL));
  const long max_rand                                = 1000000L;
  std::function<Eigen::VectorXd()> GenerateRandomVEV = [&]()
  {
    Eigen::VectorXd res(dim);
    for (auto &m : res)
      m = 1. + 299. * (rand() % max_rand) / max_rand;
    return res;
  };

  // Check if potential is invariant under group element
  std::function<int(Eigen::MatrixXd)> CheckGroupElement =
      [&](Eigen::MatrixXd GroupElement)
  {
    int result = 1;
    for (std::size_t it = 0; it < NumberOfGroupElementsChecks; it++)
    {
      auto RandomVEV = GenerateRandomVEV();
      result *= (abs(V(RandomVEV) / V(GroupElement * RandomVEV) - 1) <
                 GroupElementslMaximumRelativeError);
    }
    return result;
  };

  // Store available possible
  std::vector<Eigen::MatrixXd> C2GroupElements, SnGroupElements,
      StoreGroupElements;

  // n-plet of C2
  for (std::size_t C2Nplet = 0; C2Nplet <= dim; C2Nplet++)
  {
    // Construct bitmask. bx11...000
    std::string bitmask(C2Nplet, 1); // C2Nplet leading 1's
    bitmask.resize(dim, 0);          // NVEV - C2Nplet trailing 0's
    do                               // Permute bitmask
    {
      // Start with identity matrix
      Eigen::MatrixXd C2_element = MatrixXd::Identity(dim, dim);
      for (std::size_t i = 0; i < dim; ++i) // Flip i element
      {
        if (bitmask[i]) C2_element(i, i) = -1;
      }
      C2GroupElements.push_back(C2_element);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  }

  /*
    Currently permutations of VEV is not working

    Eigen::PermutationMatrix<Dynamic, Dynamic> PermutationMatrix(dim);
    PermutationMatrix.setIdentity();
    do
    {
      SnGroupElements.push_back(Eigen::MatrixXd(PermutationMatrix));
    } while (std::next_permutation(PermutationMatrix.indices().begin(),
                                   PermutationMatrix.indices().end()));
                                   */
  SnGroupElements.push_back(MatrixXd::Identity(dim, dim));

  // Save all group elements that keep V invariant.
  for (auto SnGroupElement : SnGroupElements)
    for (auto C2GroupElement : C2GroupElements)
    {
      if (CheckGroupElement(SnGroupElement * C2GroupElement))
        StoreGroupElements.push_back(SnGroupElement * C2GroupElement);
    }

  // Save all group elements
  this->GroupElements = StoreGroupElements;
}

std::vector<double>
MinimumTracer::ConvertToVEVDim(const std::vector<double> &point)
{
  std::vector<double> point_out;
  auto VevOrder = this->modelPointer->Get_VevOrder();

  for (std::size_t i : VevOrder)
  {
    point_out.push_back(point.at(i));
  }

  return point_out;
}

std::vector<double>
MinimumTracer::GetGlobalMinimum(const double &Temp,
                                std::vector<double> &check,
                                const std::vector<double> &start)
{
  return this->modelPointer->MinimizeOrderVEV(
      Minimizer::Minimize_gen_all(this->modelPointer,
                                  Temp,
                                  check,
                                  start,
                                  this->WhichMinimizer,
                                  this->UseMultithreading));
}

std::vector<double>
MinimumTracer::GetGlobalMinimum(const double &Temp,
                                const std::vector<double> &start)
{
  std::vector<double> check;
  return this->GetGlobalMinimum(Temp, check, start);
}

std::vector<double> MinimumTracer::GetGlobalMinimum(const double &Temp)
{
  return this->GetGlobalMinimum(
      Temp, std::vector<double>(modelPointer->get_NHiggs(), 0));
}

void MinimumTracer::IsGlobMin(Minimum &min)
{
  double num_error = 1;
  double diff      = 0;
  std::vector<double> diff_vec;

  auto glob_min =
      this->ConvertToVEVDim(this->GetGlobalMinimum(min.temp, min.point));

  ReduceVEV(glob_min);
  ConvertToNonFlatDirections(glob_min);

  auto test_min = min.point;
  ReduceVEV(test_min);
  ConvertToNonFlatDirections(test_min);

  for (std::size_t i = 0; i < glob_min.size(); i++)
  {
    diff_vec.push_back(glob_min.at(i) - test_min.at(i));
  }
  diff = L2NormVector(diff_vec);

  if (diff < num_error) // min coincides with global minimum
  {
    min.is_glob_min = true;
  }
  else // min does not coincide with global minimum
  {
    min.is_glob_min = false;
  }

  return;
}

StatusNLOStability MinimumTracer::GetStatusNLOVEV(const bool &out)
{
  std::string status;
  if (out)
  {
    return StatusNLOStability::Success;
  };
  return StatusNLOStability::NoNLOStability;
}

StatusEWSR MinimumTracer::GetStatusEWSR(const int &out)
{
  switch (out)
  {
  case 3:
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "There is EW symmetry restoration.");
    return StatusEWSR::EWSymRes;
    break;
  case 2:
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "There is no EW symmetry restoration.");
    return StatusEWSR::EWSymNonRes;
    break;
  case 1:
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "There is a flat region at high temperature. More orders are "
                  "needed to lift degeneracy.");
    return StatusEWSR::FlatRegion;
    break;
  case 0:
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "High temperature symmetry restoration calculation not "
                  "found. Check is inconclusive.");
    return StatusEWSR::Failure;
    break;
  case -1:
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Potential not bounded from below at high temperature.");
    return StatusEWSR::NotBFB;
    break;
  }
  return StatusEWSR::Failure;
}

int MinimumTracer::IsThereEWSymmetryRestoration()
{
  double T;
  double ActualSmallestEigenvalue    = 0;
  double OldSmallestEigenvalue       = 1e100;
  double EvenOlderSmallestEigenvalue = 1e200;
  double eps                         = 0.1;
  double treshold                    = 1e-6;
  double Tmax                        = 1e10;
  std::vector<double> gradient, stationary_point;
  size_t dim = this->modelPointer->get_nVEV();
  std::vector<double> point(dim, 0);

  Eigen::VectorXd GradientEigen;
  Eigen::MatrixXd HessianEigen(dim, dim);

  std::function<std::vector<double>(std::vector<double>)> dV;
  std::function<double(std::vector<double>)> V;
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian;

  Logger::Write(LoggingLevel::MinTracerDetailed,
                "Starting symmetry restoration check");

  for (double exponentT = 0; exponentT <= log(Tmax);
       exponentT += log(Tmax) / (20 * log(Tmax)))
  {
    T = exp(exponentT);
    // wrappers for potential, first and second numerical derivative
    V = [&](std::vector<double> vev)
    {
      std::vector<double> res = this->modelPointer->MinimizeOrderVEV(vev);
      if (C_UseParwani)
        return this->modelPointer->VEff(res, T) / (1 + T * T * log(T * T));
      return this->modelPointer->VEff(res, T) / (1 + T * T);
    };
    dV      = [=](auto const &arg) { return NablaNumerical(arg, V, eps); };
    Hessian = [=](auto const &arg) { return HessianNumerical(arg, V, eps); };

    ActualSmallestEigenvalue = SmallestEigenvalue(point, Hessian);

    if (abs(ActualSmallestEigenvalue / OldSmallestEigenvalue - 1) < treshold and
        abs(ActualSmallestEigenvalue / EvenOlderSmallestEigenvalue - 1) <
            treshold) // V/T^2 constant enough in T
    {
      // Save into Eigen objects
      gradient      = dV(point);
      GradientEigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
          gradient.data(), gradient.size());
      for (std::size_t i = 0; i < dim; i++)
      {
        HessianEigen.col(i) =
            Eigen::Map<Eigen::VectorXd>(Hessian(point)[i].data(), dim);
      }
      break;
    }
    EvenOlderSmallestEigenvalue = OldSmallestEigenvalue;
    OldSmallestEigenvalue       = ActualSmallestEigenvalue;
  }

  if (GradientEigen.size() == 0) return 0; // Convergence was never met

  if (ActualSmallestEigenvalue < 0) // Potential is not positively definite.
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Smallest eigenvalue at high temperature is\t" +
                      std::to_string(ActualSmallestEigenvalue));
    return -1;
  }

  if (ActualSmallestEigenvalue == 0) // Potential is not positively definite.
    return 1;

  Eigen::VectorXd StationaryPoint =
      HessianEigen.colPivHouseholderQr().solve(GradientEigen);

  // Locate global stationary point
  point = std::vector<double>(StationaryPoint.data(),
                              StationaryPoint.data() + StationaryPoint.size());

  // Check if point is reasonable.
  for (const auto &vev : point)
    if (isnan(vev)) return 0;

  // Save global minimum to use as seed point for the high temperature VEV
  HighTemperatureVEV = point;

  auto glob_min_vev_dim = modelPointer->MinimizeOrderVEV(point);
  auto EWVEV            = this->modelPointer->EWSBVEV(glob_min_vev_dim);
  if (EWVEV <= 0.5) // Symmetry restoration
  {
    return 3; // temperature for EW SR
  }
  else // (EWVEV > 0.5) EW is broken at high temperature
  {
    return 2; // SNR
  }
}

CoexPhases::CoexPhases()
{
}

CoexPhases::CoexPhases(const int &pair_id_in,
                       const Phase &false_phase_in,
                       const Phase &true_phase_in,
                       const double &Tlow_in,
                       const double &Thigh_in)
{
  coex_pair_id = pair_id_in;
  false_phase  = false_phase_in;
  true_phase   = true_phase_in;
  T_high       = Thigh_in;
  T_low        = Tlow_in;

  CalculateTc();
}

void CoexPhases::CalculateTc()
{
  Logger::Write(LoggingLevel::MinTracerDetailed,
                "Calculating critical temperature between false phase " +
                    std::to_string(false_phase.id) + " and true phase " +
                    std::to_string(true_phase.id));

  // deltaV has to be negative for the transition to occur
  auto deltaV = [&](double T)
  { return true_phase.Get(T).potential - false_phase.Get(T).potential; };

  /// Ploting deltaV
  std::stringstream ss;
  AsciiPlotter plotter("dV(T) = V(True Vacuum) - V(False Vacuum) | Phase " +
                           std::to_string(false_phase.id) + " -> Phase " +
                           std::to_string(true_phase.id),
                       100,
                       35);
  std::vector<double> plotT, plotDeltaV, plot0;
  for (double T = T_low; T <= T_high; T += (T_high - T_low) / 100)
  {
    plotT.push_back(T);
    plotDeltaV.push_back(deltaV(T));
    plot0.push_back(0);
  }
  plotter.addPlot(plotT, plot0, "", '.');
  plotter.addPlot(plotT, plotDeltaV, "dV", '*');

  plotter.xlabel("T (GeV)");
  plotter.ylabel("dV (GeV)");
  plotter.show(ss);
  Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());

  if (deltaV(T_high) > 0 and deltaV(T_low) > 0)
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "True vacuum candidate is never energetically viable.");

    crit_status = BSMPT::StatusCrit::FalseLower;
    crit_temp   = -1;
  }
  else if (deltaV(T_high) < 0 and deltaV(T_low) < 0)
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "True vacuum candidate is always energetically "
                  "viable.\nCritical temperature identified at Tc = " +
                      std::to_string(T_high) + " GeV");

    crit_status = BSMPT::StatusCrit::TrueLower;
    crit_temp   = T_high;
  }
  else if (deltaV(T_high) > 0 and deltaV(T_low) < 0)
  {
    // Binary search
    double binaryLowT    = T_low;
    double binaryMiddleT = T_low;
    double binaryHighT   = T_high;
    while (abs(binaryHighT / binaryLowT - 1) > 1e-10)
    {
      binaryMiddleT = (binaryHighT + binaryLowT) / 2;

      if (deltaV(binaryMiddleT) > 0)
      {
        binaryHighT = binaryMiddleT;
      }
      else
      {
        binaryLowT = binaryMiddleT;
      }
    }
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Critical temperature identified at Tc = " +
                      std::to_string(binaryMiddleT) + " GeV");

    crit_status = BSMPT::StatusCrit::Success;
    crit_temp   = binaryMiddleT;
  }
  else
  {
    Logger::Write(
        LoggingLevel::MinTracerDetailed,
        "False vacuum phase ends as the global minimum. New phase scan is "
        "needed!");
    crit_status = BSMPT::StatusCrit::Failure;
    crit_temp   = -2;
  }
  return;
}

std::vector<std::vector<double>>
Create1DimGrid(const std::vector<double> &point,
               const int k,
               const double low_value,
               const double high_value,
               const int nsteps)
{
  std::vector<std::vector<double>> res_vec;
  std::vector<double> res = point;
  double distance         = high_value - low_value;

  if (nsteps == 0)
  {
    res_vec.push_back(res);
    return res_vec;
  }
  else
  {
    double stepsize = distance / nsteps;

    for (int i = 0; i <= nsteps; i++)
    {
      res.at(k) = low_value + i * stepsize;
      res_vec.push_back(res);
    }
    return res_vec;
  }
}

std::vector<std::vector<double>>
Create1DimGrid(const std::vector<double> &min_start,
               const std::vector<double> &min_end,
               const int npoints)
{
  std::vector<std::vector<double>> res_vec;
  res_vec.push_back(min_start);
  auto diff_vec = min_end - min_start;
  if (npoints > 0)
  {
    auto step_vec = 1. / npoints * diff_vec;
    for (int i = 1; i <= npoints; i++)
    {
      res_vec.push_back(min_start + i * step_vec);
    }
  }
  return res_vec;
}

bool almost_the_same(const double &a,
                     const double &b,
                     const double &rel_precision,
                     const double &num_zero)
{
  if (std::abs(a) < num_zero and std::abs(b) < num_zero)
  {
    return true;
  }
  return std::abs(a - b) < std::abs(a + b) / 2 * rel_precision;
}

bool almost_the_same(const std::vector<double> &a,
                     const std::vector<double> &b,
                     const bool &allow_for_sign_flip,
                     const double &rel_precision,
                     const double &num_zero)
{
  if (a.size() != b.size())
  {
    throw std::runtime_error("Error. Vectors must have the same size.");
  }
  int count_true = 0;
  for (std::size_t i = 0; i < a.size(); i++)
  {
    if (allow_for_sign_flip)
    {
      count_true +=
          int(almost_the_same(a.at(i), b.at(i), rel_precision, num_zero));
    }
    else
    {
      count_true += int(almost_the_same(
          std::abs(a.at(i)), std::abs(b.at(i)), rel_precision, num_zero));
    }
  }
  if (std::size_t(count_true) == a.size())
  {
    return true;
  }
  else
  {
    return false;
  }
}

Phase::Phase()
{
}

Phase::Phase(const std::vector<double> &phase_start,
             const double &initialT,
             const double &finalT,
             double &globMinEndT,
             std::shared_ptr<MinimumTracer> &MinTracerIn)
{
  MinTracer = MinTracerIn;
  // Save point
  std::vector<double> phase = phase_start;
  // Reduce the VEV into the same sector
  MinTracer->ReduceVEV(phase);
  // Remove flat directions
  MinTracer->ConvertToNonFlatDirections(phase);
  std::vector<Minimum> MinimumList =
      MinTracer->TrackPhase(globMinEndT, phase, initialT, finalT);

  if (MinimumList.size() == 0) return; // Minimum tracker failed

  for (auto Min : MinimumList)
  {
    Add(Min);
  }

  T_low  = MinimumPhaseVector.front().temp;
  T_high = MinimumPhaseVector.back().temp;
}

Phase::Phase(const double &initialT,
             const std::vector<double> &phase_start,
             const double &LowT,
             const double &HighT,
             double &globMinEndT,
             std::shared_ptr<MinimumTracer> &MinTracerIn)
{
  MinTracer = MinTracerIn;
  // Save point
  std::vector<double> phase = phase_start;
  // Reduce the VEV into the same sector
  MinTracer->ReduceVEV(phase);
  // Remove flat directions
  MinTracer->ConvertToNonFlatDirections(phase);

  std::vector<Minimum> MinimumListHigh =
      MinTracer->TrackPhase(globMinEndT, phase, initialT, HighT);
  std::vector<Minimum> MinimumListLow =
      MinTracer->TrackPhase(globMinEndT, phase, initialT, LowT);

  for (auto Min : MinimumListHigh)
  {
    Add(Min);
  }
  for (auto Min : MinimumListLow)
  {
    Add(Min);
  }

  if (MinimumPhaseVector.size() == 0) return; // Found no minimum

  T_low  = MinimumPhaseVector.front().temp;
  T_high = MinimumPhaseVector.back().temp;
}

Phase::Phase(const std::vector<double> &phase_start,
             const double &initialT,
             const double &finalT,
             std::shared_ptr<MinimumTracer> &MinTracerIn)
{
  MinTracer = MinTracerIn;
  // Save point
  std::vector<double> phase = phase_start;
  // Reduce the VEV into the same sector
  MinTracer->ReduceVEV(phase);
  // Remove flat directions
  MinTracer->ConvertToNonFlatDirections(phase);

  std::vector<Minimum> MinimumList =
      MinTracer->TrackPhase(phase, initialT, finalT);

  if (MinimumList.size() == 0) return; // Minimum tracker failed

  for (auto Min : MinimumList)
  {
    Add(Min);
  }

  T_low  = MinimumPhaseVector.front().temp;
  T_high = MinimumPhaseVector.back().temp;
}

Phase::Phase(const double &initialT,
             const std::vector<double> &phase_start,
             const double &LowT,
             const double &HighT,
             std::shared_ptr<MinimumTracer> &MinTracerIn)
{
  MinTracer = MinTracerIn;
  // Save point
  std::vector<double> phase = phase_start;
  // Reduce the VEV into the same sector
  MinTracer->ReduceVEV(phase);
  // Remove flat directions
  MinTracer->ConvertToNonFlatDirections(phase);

  std::vector<Minimum> MinimumListHigh =
      MinTracer->TrackPhase(phase, initialT, HighT);
  std::vector<Minimum> MinimumListLow =
      MinTracer->TrackPhase(phase, initialT, LowT);

  for (auto Min : MinimumListHigh)
  {
    Add(Min);
  }
  for (auto Min : MinimumListLow)
  {
    Add(Min);
  }

  if (MinimumPhaseVector.size() == 0) return; // Found no minimum

  T_low  = MinimumPhaseVector.front().temp;
  T_high = MinimumPhaseVector.back().temp;
}

Phase::Phase(const double &initialT,
             const double &LowT,
             const double &HighT,
             std::shared_ptr<MinimumTracer> &MinTracerIn)
{
  MinTracer = MinTracerIn;
  std::vector<double> phase_start =
      MinTracer->ConvertToVEVDim(MinTracer->GetGlobalMinimum(initialT));

  // Reduce the VEV into the same sector
  MinTracer->ReduceVEV(phase_start);
  // Remove flat directions
  MinTracer->ConvertToNonFlatDirections(phase_start);

  if (initialT == LowT)
  {
    std::vector<Minimum> MinimumList =
        MinTracer->TrackPhase(phase_start, initialT, HighT);

    if (MinimumList.size() == 0) return; // Minimum tracker failed

    for (auto Min : MinimumList)
    {
      Add(Min);
    }
  }
  else if (initialT == HighT)
  {
    std::vector<Minimum> MinimumList =
        MinTracer->TrackPhase(phase_start, initialT, LowT);

    if (MinimumList.size() == 0) return; // Minimum tracker failed

    for (auto Min : MinimumList)
    {
      Add(Min);
    }
  }
  else if ((initialT < HighT and initialT > LowT) or
           (initialT < LowT and initialT > HighT))
  {
    std::vector<Minimum> MinimumListHigh =
        MinTracer->TrackPhase(phase_start, initialT, HighT);
    std::vector<Minimum> MinimumListLow =
        MinTracer->TrackPhase(phase_start, initialT, LowT);

    for (auto Min : MinimumListHigh)
    {
      Add(Min);
    }
    for (auto Min : MinimumListLow)
    {
      Add(Min);
    }
  }
  else
  {
    throw std::invalid_argument("Initial temperature out of bounds.");
  }

  if (MinimumPhaseVector.size() == 0) return; // Found no minimum

  T_low  = MinimumPhaseVector.front().temp;
  T_high = MinimumPhaseVector.back().temp;
}

Minimum Phase::Get(double T)
{
  if (T > T_high or T < T_low)
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Tried to Get() the minimum outside temperature range of "
                  "the phase.");
  }

  double error = 1e100;
  Minimum bestGuess;
  // Compute minimum closes to the desired temperature
  for (auto min = MinimumPhaseVector.begin(); min < MinimumPhaseVector.end();
       ++min)
  {
    if (abs(T - min->temp) < error)
    {
      error     = abs(T - min->temp);
      bestGuess = *min;
    }
  }
  if (error == 0)
  {
    // Minimum already in the list
    return bestGuess;
  }
  std::vector<Minimum> MinimumList = MinTracer->TrackPhase(
      bestGuess.point, bestGuess.temp, T, (T - bestGuess.temp), false, true);
  // Check if the TrackPhase fails
  if (MinimumList.size() > 0)
  {
    // Add new knots to the list so that Phase keeps track of the
    // new minmum information
    for (auto Min : MinimumList)
    {
      Add(Min);
    }
    // Check if last minimum found is what we wanted
    if (MinimumList.back().temp != T)
      BSMPT::Logger::Write(
          BSMPT::LoggingLevel::MinTracerDetailed,
          "Minimum tracker did not find the minimum at temperature " +
              std::to_string(T) +
              ", returning instead minimum at temperature " +
              std::to_string(MinimumList.back().temp) +
              " | dT = " + std::to_string(MinimumList.back().temp - T));
    // Return desired point
    return MinimumList.back();
  }
  return bestGuess;
}

void Phase::Add(Minimum min)
{
  // Check if phase is already there
  for (std::size_t i = 0; i < MinimumPhaseVector.size(); i++)
  {
    if (min.temp == MinimumPhaseVector[i].temp)
    {
      return;
    }
  }
  // If the list is empty add that value in.
  if (MinimumPhaseVector.size() == 0)
  {
    MinimumPhaseVector = {min};
    return;
  }
  // Temperature is the lowest yet.
  if (min.temp < MinimumPhaseVector[0].temp)
  {
    // update EdgeOfPhase
    MinimumPhaseVector.front().EdgeOfPhase = 0;
    min.EdgeOfPhase                        = -1;
    // insert min to begin of phase vector
    MinimumPhaseVector.insert(MinimumPhaseVector.begin(), min);
    T_low = min.temp;
    return;
  }
  // Is the temperature in the middle of other temperatures?
  for (std::size_t i = 0; i < MinimumPhaseVector.size() - 1; i++)
  {
    if (min.temp > MinimumPhaseVector[i].temp and
        min.temp < MinimumPhaseVector[i + 1].temp)
    {
      MinimumPhaseVector.insert(MinimumPhaseVector.begin() + i + 1, min);
      return;
    }
  }
  // If not, it is the biggest temperature yet, add it to the end and update
  // EdgeOfPhase
  MinimumPhaseVector.back().EdgeOfPhase = 0;
  min.EdgeOfPhase                       = 1;

  MinimumPhaseVector.push_back(min);
  T_high = min.temp;

  return;
}

void Vacuum::MultiStepPTMode0(const std::vector<double> &LowTempPoint_in,
                              const std::vector<double> &HighTempPoint_in)
{
  auto LowTempPoint  = LowTempPoint_in;
  auto HighTempPoint = HighTempPoint_in;

  Logger::Write(LoggingLevel::MinTracerDetailed,
                "Low-temperature phase starts at " +
                    vec_to_string(LowTempPoint) + "\n");
  Phase LowTempPhase(LowTempPoint, T_low, T_high, MinTracer);
  T_high_lowTempPhase = LowTempPhase.T_high;
  Logger::Write(LoggingLevel::MinTracerDetailed,
                "Low-temperature phase exists until T = " +
                    std::to_string(T_high_lowTempPhase) + " GeV");
  if (LowTempPhase.MinimumPhaseVector.size() > 1) addPhase(LowTempPhase);

  Logger::Write(LoggingLevel::MinTracerDetailed,
                "High-temperature phase starts at T = " +
                    vec_to_string(HighTempPoint) + " GeV");
  Phase HighTempPhase(HighTempPoint, T_high, T_low, MinTracer);
  T_low_highTempPhase = HighTempPhase.T_low;
  Logger::Write(LoggingLevel::MinTracerDetailed,
                "High-temperature phase exists until T = " +
                    std::to_string(T_low_highTempPhase) + " GeV");
  if (HighTempPhase.MinimumPhaseVector.size() > 1) addPhase(HighTempPhase);

  if (LowTempPhase.MinimumPhaseVector.size() > 1 and
      HighTempPhase.MinimumPhaseVector.size() > 1)
  {
    bool whole_temp_region_traced = T_low_highTempPhase < T_high_lowTempPhase;

    if (whole_temp_region_traced)
    {
      bool glob_min_overlap = DoGlobMinOverlap(LowTempPhase, HighTempPhase);

      if (glob_min_overlap)
      {
        Logger::Write(
            LoggingLevel::MinTracerDetailed,
            "High- and low-temperature phase coexist in temperature range [" +
                std::to_string(T_low_highTempPhase) + ", " +
                std::to_string(T_high_lowTempPhase) +
                "].\nand are "
                "still global minimum outside overlap. One-step "
                "PT identified. One-step PT mode 0 successful.");
      }
      else
      {
        Logger::Write(
            LoggingLevel::MinTracerDetailed,
            "High- and low-temperature phase coexist in temperature range [" +
                std::to_string(T_low_highTempPhase) + ", " +
                std::to_string(T_high_lowTempPhase) +
                "].\nbut are "
                "no longer global minimum outside overlap. "
                "Tracing misses the global minimum in some areas "
                "of temperature space. Try to re-run with "
                "multi-step PT mode 2 enabled.");
        status_vacuum = StatusTracing::NoGlobMinCoverage;
      }
    }
    else
    {
      Logger::Write(
          LoggingLevel::MinTracerDetailed,
          "High- and low-temperature phase do not coexist in "
          "temperature range. Possible multi-step PT detected. Try to re-run "
          "with multi-step PT mode 1 enabled.");
      status_vacuum = StatusTracing::NoCoverage;
    }
  }
  else
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "No traceable global-minimum phases found at the outer "
                  "boundaries. Abort.");
    status_vacuum = StatusTracing::NoMinsAtBoundaries;
  }
  return;
}

bool Vacuum::DoPhasesOverlap(Phase &new_phase, Phase &old_phase)
{
  bool whole_region_covered = false;

  if (new_phase.T_high >= old_phase.T_low)
  {
    if (new_phase.T_low <= old_phase.T_high)
    {
      whole_region_covered = true;
    }
  }

  return whole_region_covered;
}

bool Vacuum::DoGlobMinOverlap(const Phase &new_phase, const Phase &old_phase)
{
  bool global_minimum_overlap = false;

  if (new_phase.MinimumPhaseVector.size() > 1 and
      old_phase.MinimumPhaseVector.size() > 1)
  {
    if (new_phase.T_high >= old_phase.T_low and
        new_phase.T_low <= old_phase.T_high) // coverage check
    {
      Phase phase_low, phase_high;
      if (new_phase.T_high >= old_phase.T_high)
      {
        phase_high = new_phase;
        phase_low  = old_phase;
      }
      else if (new_phase.T_low <= old_phase.T_low)
      {
        phase_high = old_phase;
        phase_low  = new_phase;
      }

      // check for global minimum at innermost overlap
      auto min_coex_high = phase_high.Get(
          phase_low.T_high); // high-temperature phase minimum when
                             // low-temperature phase stops existing
      auto min_coex_low = phase_low.Get(
          phase_high.T_low); // low-temperature phase minimum when
                             // high-temperature phase stops existing
      MinTracer->IsGlobMin(min_coex_high);
      MinTracer->IsGlobMin(min_coex_low);
      global_minimum_overlap =
          min_coex_high.is_glob_min and min_coex_low.is_glob_min;
    }
  }
  else
  {
    global_minimum_overlap = true; // phases have no length
  }

  return global_minimum_overlap;
}

void Vacuum::MultiStepPTMode1(const std::vector<double> &LowTempPoint_in,
                              const std::vector<double> &HighTempPoint_in)
{
  auto LowTempPoint  = LowTempPoint_in;
  auto HighTempPoint = HighTempPoint_in;

  Logger::Write(LoggingLevel::MinTracerDetailed,
                "Low-temperature phase starts at T = " +
                    vec_to_string(LowTempPoint) + " GeV");
  Phase LowTempPhase(LowTempPoint, T_low, T_high, MinTracer);
  T_high_lowTempPhase = LowTempPhase.T_high;
  if (LowTempPhase.MinimumPhaseVector.size() > 1)
  {
    addPhase(LowTempPhase);
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Low-temperature phase exists until T = " +
                      std::to_string(T_high_lowTempPhase) + " GeV");
  }

  Logger::Write(LoggingLevel::MinTracerDetailed,
                "High-temperature phase starts at " +
                    vec_to_string(HighTempPoint) + "\n");
  Phase HighTempPhase(HighTempPoint, T_high, T_low, MinTracer);
  T_low_highTempPhase = HighTempPhase.T_low;
  if (HighTempPhase.MinimumPhaseVector.size() > 1)
  {
    addPhase(HighTempPhase);
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "High-temperature phase exists until T = " +
                      std::to_string(T_low_highTempPhase) + " GeV");
  }

  bool whole_temp_region_traced = DoPhasesOverlap(LowTempPhase, HighTempPhase);
  bool sides_traced             = whole_temp_region_traced;

  if (whole_temp_region_traced)
  {
    if (DoGlobMinOverlap(LowTempPhase, HighTempPhase))
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "High- and low-temperature phase coexist in "
                    "temperature range [" +
                        std::to_string(T_low_highTempPhase) + ", " +
                        std::to_string(T_high_lowTempPhase) +
                        "].\nand are "
                        "still global minimum outside overlap. One-step "
                        "PT identified. Multi-step PT mode 1 successful.");
    }
    else
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "High- and low-temperature phase overlap but are "
                    "no longer global minimum outside overlap. "
                    "Tracing misses the global minimum in some areas "
                    "of temperature space. Try to re-run with "
                    "multi-step PT mode 2 enabled.");
      status_vacuum = StatusTracing::NoGlobMinCoverage;
    }
  }
  else // check for phases in between
  {
    double tmp_low_highTempPhase = T_low_highTempPhase - 1,
           tmp_high_lowTempPhase = T_high_lowTempPhase + 1;

    auto tmp_LowTempPhase  = LowTempPhase;  // previous low temperature phase
    auto tmp_HighTempPhase = HighTempPhase; // previous high temperature phase

    bool global_minimum_overlap = true;
    bool side_global_overlap    = true;
    bool error_detected         = false;

    while (not whole_temp_region_traced and global_minimum_overlap and
           not error_detected)
    {
      Phase LowTempPhaseMiddle(tmp_high_lowTempPhase, T_low, T_high, MinTracer);
      Phase HighTempPhaseMiddle(
          tmp_low_highTempPhase, T_low, T_high, MinTracer);

      if (LowTempPhaseMiddle.MinimumPhaseVector.size() > 1)
        addPhase(LowTempPhaseMiddle);
      if (HighTempPhaseMiddle.MinimumPhaseVector.size() > 1)
        addPhase(HighTempPhaseMiddle);

      Logger::Write(
          LoggingLevel::MinTracerDetailed,
          "Intermediate phase at T = \t" +
              std::to_string(tmp_high_lowTempPhase) + " GeV exists between [ " +
              std::to_string(LowTempPhaseMiddle.T_low) + " , " +
              std::to_string(LowTempPhaseMiddle.T_high) +
              " ] GeV\nIntermediate phase at T = \t" +
              std::to_string(tmp_low_highTempPhase) + " GeV exists between [ " +
              std::to_string(HighTempPhaseMiddle.T_low) + " , " +
              std::to_string(HighTempPhaseMiddle.T_high) + " ] Gev\n");

      sides_traced = DoPhasesOverlap(tmp_LowTempPhase, LowTempPhaseMiddle) and
                     DoPhasesOverlap(tmp_HighTempPhase,
                                     HighTempPhaseMiddle); // overlap at sides

      if (not sides_traced)
      {
        Logger::Write(LoggingLevel::MinTracerDetailed,
                      "Numerical instability at phase overlap detected, no "
                      "coverage achieved.");
        status_vacuum  = StatusTracing::NoCoverage;
        error_detected = true;
      }
      else
      {
        whole_temp_region_traced =
            DoPhasesOverlap(LowTempPhaseMiddle,
                            HighTempPhaseMiddle); // innermost overlap

        if (whole_temp_region_traced)
        {
          side_global_overlap =
              DoGlobMinOverlap(tmp_LowTempPhase, LowTempPhaseMiddle) and
              DoGlobMinOverlap(tmp_HighTempPhase, HighTempPhaseMiddle);

          if (side_global_overlap)
          {
            // check for global minimum at innermost overlap
            global_minimum_overlap =
                DoGlobMinOverlap(LowTempPhaseMiddle, HighTempPhaseMiddle);

            if (global_minimum_overlap) // success
            {
              Logger::Write(
                  LoggingLevel::MinTracerDetailed,
                  "Innermost high- and low-temperature phase coexist in "
                  "temperature range "
                  "[" +
                      std::to_string(tmp_low_highTempPhase) + ", " +
                      std::to_string(tmp_high_lowTempPhase) +
                      "].\nand are "
                      "still global minimum outside overlap. Multi-step "
                      "PT identified. Multi-step PT overlap mode 1 "
                      "successful.");
            }
            else // no_glob_minimum_coverage
            {
              Logger::Write(
                  LoggingLevel::MinTracerDetailed,
                  "Innermost high- and low-temperature phase overlap but are "
                  "no longer global minimum outside overlap. "
                  "Tracing misses the global minimum in some areas "
                  "of temperature space. Try to re-run with "
                  "multi-step PT mode 2 enabled.");
              status_vacuum  = StatusTracing::NoGlobMinCoverage;
              error_detected = true;
            }
          }
          else // no_glob_minimum_coverage
          {
            Logger::Write(LoggingLevel::MinTracerDetailed,
                          "Tracing misses the global minimum in some areas "
                          "of temperature space. Try to re-run with "
                          "multi-step PT mode 2 enabled.");
            status_vacuum  = StatusTracing::NoGlobMinCoverage;
            error_detected = true;
          }
        }

        tmp_low_highTempPhase = HighTempPhaseMiddle.T_low;
        tmp_high_lowTempPhase = LowTempPhaseMiddle.T_high;

        tmp_LowTempPhase  = LowTempPhaseMiddle;
        tmp_HighTempPhase = HighTempPhaseMiddle;
      }
    }
  }

  return;
}

void Vacuum::MultiStepPTMode2(const std::vector<double> &LowTempPoint_in,
                              const std::vector<double> &HighTempPoint_in)
{
  auto LowTempPoint  = LowTempPoint_in;
  auto HighTempPoint = HighTempPoint_in;

  double T_low_newglob;
  Phase LowTempPhase(LowTempPoint, T_low, T_high, T_low_newglob, MinTracer);
  T_high_lowTempPhase = LowTempPhase.T_high;
  Logger::Write(LoggingLevel::MinTracerDetailed,
                "Low-temperature phase exists in T = [" +
                    std::to_string(LowTempPhase.T_low) + " , " +
                    std::to_string(LowTempPhase.T_high) +
                    "] GeV and is the global minimum in T = [" +
                    std::to_string(T_low) + " , " +
                    std::to_string(T_low_newglob) + "] GeV\n");
  if (LowTempPhase.MinimumPhaseVector.size() > 1) addPhase(LowTempPhase);

  double T_high_newglob;
  Phase HighTempPhase(HighTempPoint, T_high, T_low, T_high_newglob, MinTracer);
  T_low_highTempPhase = HighTempPhase.T_low;
  Logger::Write(LoggingLevel::MinTracerDetailed,
                "High-temperature phase exists until T = [" +
                    std::to_string(HighTempPhase.T_low) + " , " +
                    std::to_string(HighTempPhase.T_high) +
                    "] GeV and is the global minimum in T = [" +
                    std::to_string(T_high_newglob) + " , " +
                    std::to_string(T_high) + "] GeV\n");
  if (HighTempPhase.MinimumPhaseVector.size() > 1) addPhase(HighTempPhase);

  double deltaT = 1;

  double tmp_T_low_newglob  = (LowTempPhase.T_high > T_low_newglob)
                                  ? T_low_newglob
                                  : T_low_newglob + deltaT;
  double tmp_T_high_newglob = (HighTempPhase.T_low < T_high_newglob)
                                  ? T_high_newglob
                                  : T_high_newglob - deltaT;
  double tmp_T_low_newglob_old, tmp_T_high_newglob_old;

  while (tmp_T_low_newglob <
         tmp_T_high_newglob) // glob min not found in whole temp range
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Searching for phases in range T = [ " +
                      std::to_string(tmp_T_low_newglob) + ", " +
                      std::to_string(tmp_T_high_newglob) + " ] GeV\n");

    // Low-temperature phase
    LowTempPoint = MinTracer->ConvertToVEVDim(
        MinTracer->GetGlobalMinimum(tmp_T_low_newglob));

    Logger::Write(
        LoggingLevel::MinTracerDetailed,
        "Intermediate phase at T = " + std::to_string(tmp_T_low_newglob) +
            " GeV with " + vec_to_string(LowTempPoint) + "\n");

    tmp_T_low_newglob_old = tmp_T_low_newglob;
    Phase LowTempPhaseMiddle(tmp_T_low_newglob_old,
                             LowTempPoint,
                             T_low,
                             T_high,
                             tmp_T_low_newglob,
                             MinTracer);
    if (LowTempPhaseMiddle.T_low < LowTempPhaseMiddle.T_high)
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Intermediate phase exists between T = [ " +
                        std::to_string(LowTempPhaseMiddle.T_low) + ", " +
                        std::to_string(LowTempPhaseMiddle.T_high) + " ] GeV\n");
      addPhase(LowTempPhaseMiddle);
    }
    else
    {
      if (tmp_T_low_newglob + 1 < T_high)
        tmp_T_low_newglob += 1; // try to move out of problematic region
    }

    // High-temperature phase
    HighTempPoint = MinTracer->ConvertToVEVDim(
        MinTracer->GetGlobalMinimum(tmp_T_high_newglob));

    Logger::Write(
        LoggingLevel::MinTracerDetailed,
        "Intermediate phase at T = " + std::to_string(tmp_T_high_newglob) +
            " GeV with " + vec_to_string(HighTempPoint) + "\n");

    tmp_T_high_newglob_old = tmp_T_high_newglob;
    Phase HighTempPhaseMiddle(tmp_T_high_newglob_old,
                              HighTempPoint,
                              T_low,
                              T_high,
                              tmp_T_high_newglob,
                              MinTracer);
    if (HighTempPhaseMiddle.T_low < HighTempPhaseMiddle.T_high)
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Intermediate phase exists between T = [ " +
                        std::to_string(HighTempPhaseMiddle.T_low) + ", " +
                        std::to_string(HighTempPhaseMiddle.T_high) +
                        " ] GeV\n");
      addPhase(HighTempPhaseMiddle);
    }
    else
    {
      if (tmp_T_high_newglob > T_low + 1)
        tmp_T_high_newglob -= 1; // try to move out of problematic region
    }
  }
  return;
}

void Vacuum::PrintPhasesDiagram(int size)
{
  // If logger is disabled no need to calculate this
  if (not Logger::GetLoggingLevelStatus(LoggingLevel::MinTracerDetailed))
    return;

  std::stringstream ss;
  std::vector<double> T_list;
  std::vector<std::vector<double>> PhaseListPlot;
  std::vector<AsciiPlotter> PlotList;
  char markers[] = {".*ox-cs"};

  for (std::size_t it = 0; it < modelPointer->get_nVEV(); it++)
    PlotList.push_back(
        AsciiPlotter(modelPointer->addLegendVEV()[it], size, ceil(size / 3.)));

  for (std::size_t it = 0; it < PhasesList.size(); it++)
  {
    int nBefore =
        floor(size * (PhasesList[it].T_low - T_low) / (T_high - T_low));
    int nAfter =
        floor(size * (T_high - PhasesList[it].T_high) / (T_high - T_low));
    int nPhase = size - nBefore - nAfter;

    // Calculate VEV positions
    T_list.clear();
    PhaseListPlot.clear();
    for (int i = 0; i <= nPhase; i++)
    {
      T_list.push_back(PhasesList[it].T_low +
                       (i / double(nPhase)) *
                           (PhasesList[it].T_high - PhasesList[it].T_low));
      PhaseListPlot.push_back(PhasesList[it].Get(T_list.back()).point);
    }

    // Calculate transpose
    std::vector<std::vector<double>> PhaseListPlot_transpose(
        PhaseListPlot[0].size(), std::vector<double>(PhaseListPlot.size()));
    for (size_t i = 0; i < PhaseListPlot.size(); ++i)
      for (size_t j = 0; j < PhaseListPlot[0].size(); ++j)
        PhaseListPlot_transpose[j][i] = PhaseListPlot[i][j];

    // Do the plotting
    for (size_t i = 0; i < modelPointer->get_nVEV(); i++)
      PlotList[i].addPlot(T_list,
                          PhaseListPlot_transpose[i],
                          "Phase " + std::to_string(it),
                          markers[it % sizeof(markers)]);

    /*ss << "Phase " << std::setw(3) << it << " |" << std::string(nBefore, ' ')
       << std::string(nPhase, '.') << std::string(nAfter, ' ') << "|"
       << std::setw(7) << PhasesList[it].T_low << " - " << std::setw(7)
       << PhasesList[it].T_high << "|\n";*/

    // Print also starting and ending minima
    ss << "Phase " << std::setw(3) << it << " |" << std::string(nBefore, ' ')
       << std::string(nPhase, '.') << std::string(nAfter, ' ') << "|"
       << std::setw(7) << PhasesList[it].T_low << " - " << std::setw(7)
       << PhasesList[it].T_high << "|"
       << std::setw(2 * modelPointer->get_nVEV()) << std::setprecision(3)
       << PhasesList[it].MinimumPhaseVector.front().point << "|"
       << std::setw(2 * modelPointer->get_nVEV())
       << PhasesList[it].MinimumPhaseVector.back().point << "\n";
  }

  for (std::size_t it = 0; it < modelPointer->get_nVEV(); it++)
  {
    PlotList[it].legend();
    PlotList[it].xlabel("T (GeV)");
    PlotList[it].ylabel("Phase");
    PlotList[it].show(ss);
  }

  Logger::Write(LoggingLevel::MinTracerDetailed, ss.str());
}

Vacuum::Vacuum(const double &T_lowIn,
               const double &T_highIn,
               std::shared_ptr<MinimumTracer> &MinTracerIn,
               std::shared_ptr<Class_Potential_Origin> &modelPointerIn,
               const int &UseMultiStepPTMode,
               const int &num_pointsIn,
               const bool &do_only_tracing)
{
  T_low        = T_lowIn;
  T_high       = T_highIn;
  MinTracer    = MinTracerIn;
  modelPointer = modelPointerIn;
  num_points   = num_pointsIn;

  status_vacuum =
      StatusTracing::Success; // flipped to error code if error encountered
  status_coex_pairs = StatusCoexPair::NoCoexPairs; // flipped to success if coex
                                                   // phase pairs found

  if (UseMultiStepPTMode == -1) // default
  {
    MultiStepPTTracer(T_high);
  }
  else if (UseMultiStepPTMode >= 0)
  {
    std::vector<double> start_lowmin, start_highmin;
    for (std::size_t k = 0; k < modelPointer->get_nVEV(); k++)
    {
      start_lowmin.push_back(modelPointer->get_vevTreeMin(k));
      start_highmin.push_back(0);
    }

    std::vector<double> LowTempPoint = MinTracer->ConvertToVEVDim(
        MinTracer->GetGlobalMinimum(T_low, start_lowmin));

    if (this->modelPointer->EWSBVEV(
            this->modelPointer->MinimizeOrderVEV(LowTempPoint)) > 1e10)
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Potential non-BFB at T = 0.");
      status_vacuum = StatusTracing::Failure;
      return;
    }

    std::vector<double> HighTempPoint;
    if (MinTracer->HighTemperatureVEV.size() > 0)
    {
      HighTempPoint = MinTracer->HighTemperatureVEV;
    }
    else
    {
      HighTempPoint = MinTracer->ConvertToVEVDim(
          MinTracer->GetGlobalMinimum(T_high, start_highmin));
    }

    if (UseMultiStepPTMode == 0) // single-step phase transition mode
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Running multi-step PT mode 0.");
      MultiStepPTMode0(LowTempPoint, HighTempPoint);
    }
    else if (UseMultiStepPTMode == 1) // enforce tracing coverage
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Running multi-step PT mode 1.");
      MultiStepPTMode1(LowTempPoint, HighTempPoint);
    }
    else if (UseMultiStepPTMode == 2) // enforce global minimum tracing coverage
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Running multi-step PT mode 2.");
      MultiStepPTMode2(LowTempPoint, HighTempPoint);
    }
    else if (UseMultiStepPTMode == 3) // automatic mode
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Running multi-step PT mode auto.");
      MultiStepPTMode1(LowTempPoint, HighTempPoint);
      if (status_vacuum == StatusTracing::NoGlobMinCoverage)
      {
        status_vacuum = StatusTracing::Success;
        MultiStepPTMode2(LowTempPoint, HighTempPoint);
      }
    }
  }

  // trace EW minimum (still at least a local minimum)
  Minimum EWMin;
  EWMin.temp  = 0;
  EWMin.point = modelPointer->get_vevTreeMin();
  if (MinimumFoundAlready(EWMin) == -1)
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Tracking the EW minimum at " + vec_to_string(EWMin.point));
    Phase EWphase(EWMin.point, T_low, T_high, MinTracer);
    addPhase(EWphase);
    print(EWphase);
  }
  else
    Logger::Write(LoggingLevel::MinTracerDetailed, "EW minimum already found.");

  // VEVs found from minima splitting
  std::vector<Minimum> SavedMinimaFromVEVSplitting_Copy =
      MinTracer->SavedMinimaFromVEVSplitting; // Make a copy of the found phases
  for (auto Min : SavedMinimaFromVEVSplitting_Copy)
  {
    if (MinimumFoundAlready(Min) > -1) continue;
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Found a new phase from VEV splitting at T = " +
                      std::to_string(Min.temp));
    Phase newPhase(Min.temp, Min.point, T_low, T_high, MinTracer);
    addPhase(newPhase);
    print(newPhase);
  }

  if (PhasesList.size() == 0) // no phases could be found
  {
    status_vacuum = StatusTracing::Failure;
  }

  if (status_vacuum == StatusTracing::Success or
      (UseMultiStepPTMode != 0 and
       status_vacuum == StatusTracing::NoCoverage)) // no_coverage can get fixed
                                                    // in setCoexRegion
  {
    // sort phases in decending T_high
    std::sort(PhasesList.begin(),
              PhasesList.end(),
              [](auto a, auto b) { return a.T_high > b.T_high; });

    // assign ids to phases
    for (std::size_t i = 0; i < PhasesList.size(); i++)
    {
      PhasesList[i].id = i;
    }

    // identify coexisiting phase regions
    setCoexRegion(UseMultiStepPTMode); // can flip status_vacuum to error code

    if (PhasesList.size() > 0)
    {
      PrintPhasesDiagram();
    }

    if ((status_coex_pairs == StatusCoexPair::Success) and
        (not do_only_tracing))
    {
      // identify coexisting phase pairs
      setCoexPhases();
    }
  }

  return;
}

void Vacuum::MultiStepPTTracer(const double &Temp, const double &deltaT)
{
  if (Temp <= T_low)
  {
    auto glob_min = MinTracer->GetGlobalMinimum(T_low);
    if (this->modelPointer->EWSBVEV(glob_min) > 1e10)
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Potential non-BFB at T = 0.");
      status_vacuum = StatusTracing::Failure;
    }
    else
    {
      // low temperature phase
      Phase phase(T_low, T_low, T_high, MinTracer);
      addPhase(phase);
      print(phase);

      // test equally-spaced point grid
      for (int i = 1; i <= num_points; i++)
      {
        Minimum min;
        min.temp = T_low + (T_high - T_low) / (num_points + 1) * i;
        min.point =
            MinTracer->ConvertToVEVDim(MinTracer->GetGlobalMinimum(min.temp));
        MinTracer->ReduceVEV(min.point);
        MinTracer->ConvertToNonFlatDirections(min.point);

        if (MinimumFoundAlready(min) == -1) // found new phase
        {
          Logger::Write(
              LoggingLevel::MinTracerDetailed,
              "-------------------------------------------------------");
          Phase inter_phase(min.temp, T_high, T_low, MinTracer);
          addPhase(inter_phase);
          print(inter_phase);
        }
        else
        {
          Logger::Write(LoggingLevel::MinTracerDetailed,
                        "Point at T = " + std::to_string(min.temp) +
                            " GeV belongs to already traced phase.");
        }
      }
    }
    return;
  }
  else if ((Temp == T_high) and (MinTracer->HighTemperatureVEV.size() > 0))
  {
    Minimum HTMin;
    HTMin.temp  = T_high;
    HTMin.point = MinTracer->HighTemperatureVEV;
    if (MinimumFoundAlready(HTMin) == -1)
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Tracking the high-temperature minimum at " +
                        vec_to_string(HTMin.point));
      Phase HTphase(HTMin.point, T_high, T_low, MinTracer);
      addPhase(HTphase);
      print(HTphase);
      MultiStepPTTracer(HTphase.T_low, -1);
      return;
    }
  }

  Minimum min;
  min.temp  = Temp;
  min.point = MinTracer->ConvertToVEVDim(MinTracer->GetGlobalMinimum(Temp));

  // Reduce the VEV into the same sector
  MinTracer->ReduceVEV(min.point);
  // Remove flat directions
  MinTracer->ConvertToNonFlatDirections(min.point);

  if (MinimumFoundAlready(min) == -1) // found new phase
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "-------------------------------------------------------");
    Phase phase(Temp, T_high, T_low, MinTracer);
    addPhase(phase);
    print(phase);

    MultiStepPTTracer(phase.T_low, -1);
  }
  else // found last traced phase again
  {
    // Signal of phase splitting
    MultiStepPTTracer(Temp + deltaT);
  }

  return;
}

void Vacuum::print(const Phase &phase)
{
  if (phase.MinimumPhaseVector.size() > 1)
  {
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "New phase found between " + std::to_string(phase.T_low) +
                      " GeV and " + std::to_string(phase.T_high) + " GeV");
  }
  return;
}

void Vacuum::setCoexPhases()
{
  std::stringstream ss1, ss2;
  int pair_id = 0;

  for (std::size_t i = 0; i < PhasesList.size(); i++)
  {
    auto false_phase = PhasesList.at(i);
    ss1 << "Phase " << false_phase.id << " exists between T = ["
        << false_phase.T_low << ", " << false_phase.T_high << "] GeV\n";
    for (std::size_t j = i + 1; j < PhasesList.size(); j++)
    {
      auto true_phase = PhasesList.at(j);
      if (true_phase.T_high >= false_phase.T_low and
          true_phase.T_low <= false_phase.T_high)
      {
        double T_low_overlap  = max(false_phase.T_low, true_phase.T_low);
        double T_high_overlap = min(false_phase.T_high, true_phase.T_high);

        CoexPhases phase_pair(
            pair_id, false_phase, true_phase, T_low_overlap, T_high_overlap);
        CoexPhasesList.push_back(phase_pair);

        ss2 << "Pair " << pair_id << " (" << false_phase.id << ", "
            << true_phase.id << ") with T = [" << T_low_overlap << ", "
            << T_high_overlap << "] GeV and Tc = " << phase_pair.crit_temp
            << " (" << StatusCritToString.at(phase_pair.crit_status) << ")\n";

        pair_id += 1;
      }
    }
  }

  Logger::Write(LoggingLevel::MinTracerDetailed, ss1.str());
  Logger::Write(LoggingLevel::MinTracerDetailed, ss2.str());
}

void Vacuum::setCoexRegion(const int &UseMultiStepPTMode)
{
  std::vector<Minimum> edgesList, edgesListResult;
  std::vector<double> tempList;
  Logger::Write(LoggingLevel::MinTracerDetailed,
                "Total number of phases identified: " +
                    std::to_string(PhasesList.size()));

  if (PhasesList.size() > 0)
  {
    for (auto i : PhasesList)
    {
      edgesList.push_back(i.MinimumPhaseVector.front());
      edgesList.push_back(i.MinimumPhaseVector.back());
    }

    for (auto i : edgesList)
    {
      tempList.push_back(i.temp);
    }

    std::sort(
        tempList.begin(), tempList.end(), [](auto a, auto b) { return a > b; });
    tempList.erase(unique(tempList.begin(), tempList.end()), tempList.end());

    int EdgeOfPhaseatTemp = 0;
    for (auto temp : tempList)
    {
      for (auto edge : edgesList)
      {
        if (edge.temp == temp)
        {
          EdgeOfPhaseatTemp += edge.EdgeOfPhase;
        }
      }
      Minimum min;
      min.temp        = temp;
      min.EdgeOfPhase = EdgeOfPhaseatTemp;
      edgesListResult.push_back(min);
    }

    // order list decending in temperature
    std::sort(edgesList.begin(),
              edgesList.end(),
              [](auto a, auto b) { return a.temp > b.temp; });
    std::sort(edgesListResult.begin(),
              edgesListResult.end(),
              [](auto a, auto b) { return a.temp > b.temp; });

    int numPhases = 0;

    bool no_gap_found = true;

    for (std::size_t i = 0; i < edgesListResult.size() - 1; i++)
    {
      numPhases = edgesListResult[i].EdgeOfPhase;
      if (numPhases > 1)
      {
        status_coex_pairs = StatusCoexPair::Success;
      }
      else if (numPhases <= 0) // found a non-traced gap in temperature
      {
        no_gap_found       = false;
        double T_low_hole  = edgesListResult[i + 1].temp;
        double T_high_hole = edgesListResult[i].temp;

        if (T_low_hole + 1e-6 < T_high_hole) // threshold set to 1e-6 GeV
        {
          Logger::Write(LoggingLevel::MinTracerDetailed,
                        "\nThere are phases missing between " +
                            std::to_string(T_low_hole) + " GeV and " +
                            std::to_string(T_high_hole) + " GeV!");
          status_vacuum = StatusTracing::NoCoverage;

          if (not(UseMultiStepPTMode == 0))
          {
            Logger::Write(LoggingLevel::MinTracerDetailed,
                          "\nTry to patch up gap between " +
                              std::to_string(T_low_hole) + " GeV and " +
                              std::to_string(T_high_hole) + " GeV.");

            // try to patch up holes in tracing
            Minimum min;
            min.temp  = (T_high_hole + T_low_hole) / 2;
            min.point = MinTracer->ConvertToVEVDim(
                MinTracer->GetGlobalMinimum(min.temp));
            MinTracer->ReduceVEV(min.point);
            MinTracer->ConvertToNonFlatDirections(min.point);
            Phase inter_phase(min.temp, T_high, T_low, MinTracer);
            addPhase(inter_phase);
            print(inter_phase);
            // If we had sucess finding a new phase.
            if (inter_phase.MinimumPhaseVector.size() >
                2) // more than just endpoints found
            {
              status_vacuum = StatusTracing::Success;
              setCoexRegion(UseMultiStepPTMode);
            }
          }
        }
        else
        {
          status_vacuum = StatusTracing::NoCoverage;
        }
      }
    }

    if (no_gap_found) // correct status code in case local EW minimum covers up
                      // range
    {
      if (status_vacuum == StatusTracing::NoCoverage)
      {
        status_vacuum = StatusTracing::NoGlobMinCoverage;
      }
    }
  }

  return;
}

void Vacuum::addPhase(Phase &phase)
{
  if (phase.MinimumPhaseVector.size() <= 1)
  {
    if (phase.MinimumPhaseVector.size() > 0)
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    vec_to_string(phase.MinimumPhaseVector.at(0).point));
    Logger::Write(LoggingLevel::MinTracerDetailed,
                  "Phase was empty. Nothing was added.");
    return;
  }
  // If the step size is too high and the T_high and T_low are too close we
  // might have a phase with only 2 minima. This is not ideal so we look for
  // more in the middle.
  if (phase.MinimumPhaseVector.size() == 2)
  {
    // If only two minima were found look for more.
    std::vector<Minimum> MinimumList = MinTracer->TrackPhase(
        phase.MinimumPhaseVector[0].point,
        phase.MinimumPhaseVector[0].temp,
        phase.MinimumPhaseVector[1].temp,
        (phase.MinimumPhaseVector[1].temp - phase.MinimumPhaseVector[0].temp) /
            3,
        false);
    // Check if the TrackPhase fails
    if (MinimumList.size() > 0)
    {
      // Add new knots to the list so that Phase keeps track of the
      // new minumum information
      for (auto Min : MinimumList)
      {
        phase.Add(Min);
      }
    }
    else
    {
      Logger::Write(LoggingLevel::MinTracerDetailed,
                    "Phase is too unstable. Throwing it away.");
      return;
    }
  }
  for (auto &existingPhase : PhasesList)
  {
    // Calculate the maximum distance of phase to each other phases
    double Temp_min = std::max(existingPhase.T_low, phase.T_low);
    double Temp_max = std::min(existingPhase.T_high, phase.T_high);

    if (Temp_min >= Temp_max)
    {
      continue; // Phases never exist at the same time
    }
    else
    {
      double dT           = (Temp_max - Temp_min);
      double max_distance = 0;
      double numSteps     = 12;
      for (double n = 1; n < numSteps; n++)
      {
        max_distance = std::max(
            max_distance,
            L2NormVector(
                existingPhase.Get(Temp_min + (n / numSteps) * dT).point -
                phase.Get(Temp_min + (n / numSteps) * dT).point));
      }
      if (max_distance < 1)
      {
        Logger::Write(
            LoggingLevel::MinTracerDetailed,
            "The phase starting at T = " + std::to_string(phase.T_low) +
                " GeV and ending at T = " + std::to_string(phase.T_high) +
                " GeV could not be added because it coincides with an "
                "already "
                "existing phase that starts at T = " +
                std::to_string(existingPhase.T_low) + " GeV and ends at T = " +
                std::to_string(existingPhase.T_high) + " GeV.");

        // Add Minimum of other phase to the already existing phase.
        for (auto min : phase.MinimumPhaseVector)
          existingPhase.Add(min);
        return; // The phases coincide. Abort!
      }
    }
  }

  // Set starting of phase
  phase.MinimumPhaseVector.front().EdgeOfPhase = -1;
  // Set ending of phase
  phase.MinimumPhaseVector.back().EdgeOfPhase = 1;
  PhasesList.push_back(phase);
  return;
}

int Vacuum::MinimumFoundAlready(const Minimum &minimum)
{
  double MaximumAllowedError = 1;
  for (std::size_t ind = 0; ind < PhasesList.size(); ind++)
  {
    // Check if it is outside temperature range
    if (minimum.temp > PhasesList[ind].T_high or
        minimum.temp < PhasesList[ind].T_low)
      continue;

    std::vector<double> Distance;
    std::transform(minimum.point.begin(),
                   minimum.point.end(),
                   PhasesList[ind].Get(minimum.temp).point.begin(),
                   std::back_inserter(Distance),
                   std::minus<double>());

    if (L2NormVector(Distance) / this->modelPointer->get_nVEV() <
        MaximumAllowedError)
      return ind;
  }
  return -1;
}

std::vector<std::string> MinimumTracer::GetLegend(const int &num_coex_phases,
                                                  const bool &do_gw_calc)
{
  std::vector<std::string> legend;

  legend.push_back("status_nlo_stability");
  legend.push_back("status_ewsr");
  legend.push_back("status_tracing");
  legend.push_back("status_coex_pairs");
  legend.push_back("runtime");

  for (int i = 0; i < num_coex_phases; i++)
  {
    legend.push_back("status_crit_" + std::to_string(i));
    legend.push_back("T_crit_" + std::to_string(i));
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_crit_false_" + std::to_string(i)));
    }
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_crit_true_" + std::to_string(i)));
    }
    legend.push_back("status_bounce_sol_" + std::to_string(i));
    legend.push_back("status_nucl_approx_" + std::to_string(i));
    legend.push_back("T_nucl_approx_" + std::to_string(i));
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_nucl_approx_false_" + std::to_string(i)));
    }
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_nucl_approx_true_" + std::to_string(i)));
    }
    legend.push_back("status_nucl_" + std::to_string(i));
    legend.push_back("T_nucl_" + std::to_string(i));
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_nucl_false_" + std::to_string(i)));
    }
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_nucl_true_" + std::to_string(i)));
    }
    legend.push_back("status_perc_" + std::to_string(i));
    legend.push_back("T_perc_" + std::to_string(i));
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_perc_false_" + std::to_string(i)));
    }
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_perc_true_" + std::to_string(i)));
    }
    legend.push_back("status_compl_" + std::to_string(i));
    legend.push_back("T_compl_" + std::to_string(i));
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_compl_false_" + std::to_string(i)));
    }
    for (std::size_t j = 0; j < this->modelPointer->addLegendVEV().size(); j++)
    {
      legend.push_back(this->modelPointer->addLegendVEV().at(j).append(
          "_compl_true_" + std::to_string(i)));
    }
    if (do_gw_calc)
    {
      legend.push_back("status_gw_" + std::to_string(i));
      legend.push_back("T_star_" + std::to_string(i));
      legend.push_back("T_reh_" + std::to_string(i));
      legend.push_back("v_wall_" + std::to_string(i));
      legend.push_back("alpha_PT_" + std::to_string(i));
      legend.push_back("beta/H_" + std::to_string(i));
      legend.push_back("kappa_col_" + std::to_string(i));
      legend.push_back("kappa_sw_" + std::to_string(i));
      legend.push_back("eps_turb_" + std::to_string(i));
      legend.push_back("cs_f_" + std::to_string(i));
      legend.push_back("cs_t_" + std::to_string(i));
      legend.push_back("fb_col_" + std::to_string(i));
      legend.push_back("h2Omegab_col_" + std::to_string(i));
      legend.push_back("f_1_sw_" + std::to_string(i));
      legend.push_back("f_2_sw_" + std::to_string(i));
      legend.push_back("h2Omega_2_sw_" + std::to_string(i));
      legend.push_back("f_1_turb_" + std::to_string(i));
      legend.push_back("f_2_turb_" + std::to_string(i));
      legend.push_back("h2Omega_2_turb_" + std::to_string(i));
      legend.push_back("SNR(LISA-3yrs)_col_" + std::to_string(i));
      legend.push_back("SNR(LISA-3yrs)_sw_" + std::to_string(i));
      legend.push_back("SNR(LISA-3yrs)_turb_" + std::to_string(i));
      legend.push_back("SNR(LISA-3yrs)_" + std::to_string(i));
    }
  }

  legend.push_back("transition_history");

  return legend;
}

} // namespace BSMPT
