// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/utility/utility.h>

TEST_CASE("Check vector . vector product", "[utility]")
{
  using namespace BSMPT;
  std::vector<double> v1 = {2.237, 0.189, -1.202, 6.911, 2.754};
  std::vector<double> v2 = {8.39, -3.968, 8.046, 3.925, 2.209};
  REQUIRE(v1 * v2 == Approx(41.556447).margin(1e-10));
}

TEST_CASE("Check matrix . vector", "[utility]")
{
  using namespace BSMPT;
  std::vector<std::vector<double>> H = {
      {-3.493, -1.695, -8.702, -3.008, -3.011},
      {7.26, 9.379, -2.769, -6.25, 7.548},
      {4.295, -1.01, -7.67, -4.673, -2.335},
      {-6.258, 4.273, 1.702, -3.612, -8.57},
      {-5.182, -5.749, 2.92, 0.864, 4.06}};
  std::vector<double> v1 = {0.856, 8.936, -5.996, 4.987, 9.301};

  std::vector<double> result = H * v1;

  REQUIRE(result[0] == Approx(-8.965543).margin(1e-10));
  REQUIRE(result[1] == Approx(145.6634260).margin(1e-10));
  REQUIRE(result[2] == Approx(-4.381606).margin(1e-10));
  REQUIRE(result[3] == Approx(-75.101126).margin(1e-10));
  REQUIRE(result[4] == Approx(-31.246348).margin(1e-10));
}

TEST_CASE("Check vector . matrix . vector", "[utility]")
{
  using namespace BSMPT;
  std::vector<std::vector<double>> H = {{-6.945, -6.422, 8.201, 7.182, -4.665},
                                        {9.531, 9.744, 2.673, -5.752, -6.406},
                                        {-3.209, 8.142, -3.635, -7.9, 6.279},
                                        {-7.554, -3.314, 5.775, 2.42, -5.703},
                                        {-3.549, 2.838, 6.843, 3.395, 2.009}};
  std::vector<double> v1             = {-8.429, -9.61, -1.958, 8.316, -7.775};
  std::vector<double> v2             = {4.882, -0.446, 9.56, -7.453, 0.934};

  REQUIRE(v1 * (H * v2) == Approx(-1149.174865593).margin(1e-10));
}