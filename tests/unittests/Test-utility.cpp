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

TEST_CASE("Check Transpose", "[utility]")
{
  using namespace BSMPT;
  std::vector<std::vector<double>> H = {
      {-7.697, 3.259, -2.48, 9.655, -2.463},
      {1.17, -1.936, -0.053, 9.129, 4.828},
      {-4.173, -1.695, -8.372, -8.434, -7.375},
      {-1.807, -3.037, 4.649, 5.487, -5.915},
      {-2.805, -8.22, 3.369, -4.824, -7.272}};

  std::vector<std::vector<double>> H_Transpose = {
      {{-7.697, 1.17, -4.173, -1.807, -2.805},
       {3.259, -1.936, -1.695, -3.037, -8.22},
       {-2.48, -0.053, -8.372, 4.649, 3.369},
       {9.655, 9.129, -8.434, 5.487, -4.824},
       {-2.463, 4.828, -7.375, -5.915, -7.272}}};

  REQUIRE(H_Transpose == BSMPT::Transpose(H));
}

TEST_CASE("Check L2NormVector", "[utility]")
{
  using namespace BSMPT;
  std::vector<double> v1 = {0.856, 8.936, -5.996, 4.987, 9.301};

  REQUIRE(BSMPT::L2NormVector(v1) == Approx(15.096874).margin(1e-10));
}

TEST_CASE("Check Li2 function", "[utility]")
{
  using namespace BSMPT;
  // https://en.wikipedia.org/wiki/Dilogarithm
  REQUIRE(Li2(-1) == Approx(-pow(M_PI, 2) / 12.).margin(1e-10));
  REQUIRE(Li2(0) == Approx(0).margin(1e-10));
  REQUIRE(Li2(1. / 2.) ==
          Approx(pow(M_PI, 2) / 12. - pow(log(2), 2) / 2).margin(1e-10));
  REQUIRE(Li2((1 - sqrt(5)) / 2) ==
          Approx(-pow(M_PI, 2) / 15. + pow(log((sqrt(5) + 1) / 2), 2) / 2.)
              .margin(1e-10));
  REQUIRE(Li2(-(1 + sqrt(5)) / 2) ==
          Approx(-pow(M_PI, 2) / 10. - pow(log((sqrt(5) + 1) / 2), 2))
              .margin(1e-10));
  REQUIRE(Li2((3 - sqrt(5)) / 2) ==
          Approx(pow(M_PI, 2) / 15. - pow(log((sqrt(5) + 1) / 2), 2))
              .margin(1e-10));
  REQUIRE(Li2((-1 + sqrt(5)) / 2) ==
          Approx(pow(M_PI, 2) / 10. - pow(log((sqrt(5) + 1) / 2), 2))
              .margin(1e-10));

  REQUIRE(Li2(1 - 1e-5) == Approx(1.644808937).margin(1e-10));
  REQUIRE(Li2(1 - 1e-4) == Approx(1.643912984).margin(1e-10));
  REQUIRE(Li2(1 - 1e-3) == Approx(1.637022605).margin(1e-10));
  REQUIRE(Li2(1 - 1e-2) == Approx(1.588625448).margin(1e-10));
  REQUIRE(Li2(1 - 1e-1) == Approx(1.299714723).margin(1e-10));
  REQUIRE(Li2(1 - 1e0) == Approx(0.).margin(1e-10));
  REQUIRE(Li2(1 - 1e1) == Approx(-3.950663778).margin(1e-10));
  REQUIRE(Li2(1 - 1e2) == Approx(-12.19242167).margin(1e-10));
  REQUIRE(Li2(1 - 1e3) == Approx(-25.4955641).margin(1e-10));
  REQUIRE(Li2(1 - 1e4) == Approx(-44.05909787).margin(1e-10));
  REQUIRE(Li2(1 - 1e5) == Approx(-67.91853532).margin(1e-10));

  REQUIRE(Li2(-1 + 1e-5) == Approx(-0.8224601019).margin(1e-10));
  REQUIRE(Li2(-1 + 1e-4) == Approx(-0.8223977177).margin(1e-10));
  REQUIRE(Li2(-1 + 1e-3) == Approx(-0.8217737896).margin(1e-10));
  REQUIRE(Li2(-1 + 1e-2) == Approx(-0.8155258815).margin(1e-10));
  REQUIRE(Li2(-1 + 1e-1) == Approx(-0.7521631792).margin(1e-10));

  REQUIRE(Li2(1e-5) == Approx(0.000010000025).margin(1e-10));
  REQUIRE(Li2(1e-4) == Approx(0.0001000025001).margin(1e-10));
  REQUIRE(Li2(1e-3) == Approx(0.001000250111).margin(1e-10));
  REQUIRE(Li2(1e-2) == Approx(0.01002511174).margin(1e-10));
  REQUIRE(Li2(1e-1) == Approx(0.1026177911).margin(1e-10));

  REQUIRE(Li2(-1e-5) == Approx(-9.999975e-6).margin(1e-10));
  REQUIRE(Li2(-1e-4) == Approx(-0.00009999750025).margin(1e-10));
  REQUIRE(Li2(-1e-3) == Approx(-0.000999750111).margin(1e-10));
  REQUIRE(Li2(-1e-2) == Approx(-0.00997511049).margin(1e-10));
  REQUIRE(Li2(-1e-1) == Approx(-0.09760523523).margin(1e-10));
  REQUIRE(Li2(-1e0) == Approx(-0.8224670334).margin(1e-10));
  REQUIRE(Li2(-1e1) == Approx(-4.198277887).margin(1e-10));
  REQUIRE(Li2(-1e2) == Approx(-12.23875518).margin(1e-10));
  REQUIRE(Li2(-1e3) == Approx(-25.50247581).margin(1e-10));
  REQUIRE(Li2(-1e4) == Approx(-44.06001895).margin(1e-10));
  REQUIRE(Li2(-1e5) == Approx(-67.91865045).margin(1e-10));

  REQUIRE(Li2(-1.) == Approx(-0.8224670334).margin(1e-10));
  REQUIRE(Li2(-0.99) == Approx(-0.8155258815).margin(1e-10));
  REQUIRE(Li2(-0.98) == Approx(-0.8085652776).margin(1e-10));
  REQUIRE(Li2(-0.97) == Approx(-0.801585083).margin(1e-10));
  REQUIRE(Li2(-0.96) == Approx(-0.7945851575).margin(1e-10));
  REQUIRE(Li2(-0.95) == Approx(-0.7875653589).margin(1e-10));
  REQUIRE(Li2(-0.94) == Approx(-0.7805255435).margin(1e-10));
  REQUIRE(Li2(-0.93) == Approx(-0.773465566).margin(1e-10));
  REQUIRE(Li2(-0.92) == Approx(-0.7663852791).margin(1e-10));
  REQUIRE(Li2(-0.91) == Approx(-0.7592845337).margin(1e-10));
  REQUIRE(Li2(-0.9) == Approx(-0.7521631792).margin(1e-10));
  REQUIRE(Li2(-0.89) == Approx(-0.7450210629).margin(1e-10));
  REQUIRE(Li2(-0.88) == Approx(-0.7378580301).margin(1e-10));
  REQUIRE(Li2(-0.87) == Approx(-0.7306739245).margin(1e-10));
  REQUIRE(Li2(-0.86) == Approx(-0.7234685878).margin(1e-10));
  REQUIRE(Li2(-0.85) == Approx(-0.7162418594).margin(1e-10));
  REQUIRE(Li2(-0.84) == Approx(-0.7089935771).margin(1e-10));
  REQUIRE(Li2(-0.83) == Approx(-0.7017235764).margin(1e-10));
  REQUIRE(Li2(-0.82) == Approx(-0.6944316907).margin(1e-10));
  REQUIRE(Li2(-0.81) == Approx(-0.6871177515).margin(1e-10));
  REQUIRE(Li2(-0.8) == Approx(-0.6797815878).margin(1e-10));
  REQUIRE(Li2(-0.79) == Approx(-0.6724230268).margin(1e-10));
  REQUIRE(Li2(-0.78) == Approx(-0.6650418931).margin(1e-10));
  REQUIRE(Li2(-0.77) == Approx(-0.6576380092).margin(1e-10));
  REQUIRE(Li2(-0.76) == Approx(-0.6502111952).margin(1e-10));
  REQUIRE(Li2(-0.75) == Approx(-0.6427612688).margin(1e-10));
  REQUIRE(Li2(-0.74) == Approx(-0.6352880455).margin(1e-10));
  REQUIRE(Li2(-0.73) == Approx(-0.6277913381).margin(1e-10));
  REQUIRE(Li2(-0.72) == Approx(-0.620270957).margin(1e-10));
  REQUIRE(Li2(-0.71) == Approx(-0.61272671).margin(1e-10));
  REQUIRE(Li2(-0.7) == Approx(-0.6051584023).margin(1e-10));
  REQUIRE(Li2(-0.69) == Approx(-0.5975658366).margin(1e-10));
  REQUIRE(Li2(-0.68) == Approx(-0.5899488126).margin(1e-10));
  REQUIRE(Li2(-0.67) == Approx(-0.5823071275).margin(1e-10));
  REQUIRE(Li2(-0.66) == Approx(-0.5746405756).margin(1e-10));
  REQUIRE(Li2(-0.65) == Approx(-0.5669489483).margin(1e-10));
  REQUIRE(Li2(-0.64) == Approx(-0.5592320342).margin(1e-10));
  REQUIRE(Li2(-0.63) == Approx(-0.5514896188).margin(1e-10));
  REQUIRE(Li2(-0.62) == Approx(-0.5437214845).margin(1e-10));
  REQUIRE(Li2(-0.61) == Approx(-0.5359274109).margin(1e-10));
  REQUIRE(Li2(-0.6) == Approx(-0.528107174).margin(1e-10));
  REQUIRE(Li2(-0.59) == Approx(-0.520260547).margin(1e-10));
  REQUIRE(Li2(-0.58) == Approx(-0.5123872996).margin(1e-10));
  REQUIRE(Li2(-0.57) == Approx(-0.504487198).margin(1e-10));
  REQUIRE(Li2(-0.56) == Approx(-0.4965600052).margin(1e-10));
  REQUIRE(Li2(-0.55) == Approx(-0.4886054807).margin(1e-10));
  REQUIRE(Li2(-0.54) == Approx(-0.4806233803).margin(1e-10));
  REQUIRE(Li2(-0.53) == Approx(-0.4726134562).margin(1e-10));
  REQUIRE(Li2(-0.52) == Approx(-0.4645754568).margin(1e-10));
  REQUIRE(Li2(-0.51) == Approx(-0.4565091268).margin(1e-10));
  REQUIRE(Li2(-0.5) == Approx(-0.4484142069).margin(1e-10));
  REQUIRE(Li2(-0.49) == Approx(-0.440290434).margin(1e-10));
  REQUIRE(Li2(-0.48) == Approx(-0.4321375408).margin(1e-10));
  REQUIRE(Li2(-0.47) == Approx(-0.4239552559).margin(1e-10));
  REQUIRE(Li2(-0.46) == Approx(-0.4157433035).margin(1e-10));
  REQUIRE(Li2(-0.45) == Approx(-0.4075014037).margin(1e-10));
  REQUIRE(Li2(-0.44) == Approx(-0.3992292721).margin(1e-10));
  REQUIRE(Li2(-0.43) == Approx(-0.3909266197).margin(1e-10));
  REQUIRE(Li2(-0.42) == Approx(-0.3825931528).margin(1e-10));
  REQUIRE(Li2(-0.41) == Approx(-0.3742285731).margin(1e-10));
  REQUIRE(Li2(-0.4) == Approx(-0.3658325775).margin(1e-10));
  REQUIRE(Li2(-0.39) == Approx(-0.3574048577).margin(1e-10));
  REQUIRE(Li2(-0.38) == Approx(-0.3489451006).margin(1e-10));
  REQUIRE(Li2(-0.37) == Approx(-0.3404529876).margin(1e-10));
  REQUIRE(Li2(-0.36) == Approx(-0.331928195).margin(1e-10));
  REQUIRE(Li2(-0.35) == Approx(-0.3233703936).margin(1e-10));
  REQUIRE(Li2(-0.34) == Approx(-0.3147792486).margin(1e-10));
  REQUIRE(Li2(-0.33) == Approx(-0.3061544195).margin(1e-10));
  REQUIRE(Li2(-0.32) == Approx(-0.2974955599).margin(1e-10));
  REQUIRE(Li2(-0.31) == Approx(-0.2888023175).margin(1e-10));
  REQUIRE(Li2(-0.3) == Approx(-0.2800743338).margin(1e-10));
  REQUIRE(Li2(-0.29) == Approx(-0.2713112439).margin(1e-10));
  REQUIRE(Li2(-0.28) == Approx(-0.2625126766).margin(1e-10));
  REQUIRE(Li2(-0.27) == Approx(-0.2536782541).margin(1e-10));
  REQUIRE(Li2(-0.26) == Approx(-0.2448075917).margin(1e-10));
  REQUIRE(Li2(-0.25) == Approx(-0.2359002977).margin(1e-10));
  REQUIRE(Li2(-0.24) == Approx(-0.2269559734).margin(1e-10));
  REQUIRE(Li2(-0.23) == Approx(-0.2179742128).margin(1e-10));
  REQUIRE(Li2(-0.22) == Approx(-0.2089546022).margin(1e-10));
  REQUIRE(Li2(-0.21) == Approx(-0.1998967202).margin(1e-10));
  REQUIRE(Li2(-0.2) == Approx(-0.1908001378).margin(1e-10));
  REQUIRE(Li2(-0.19) == Approx(-0.1816644174).margin(1e-10));
  REQUIRE(Li2(-0.18) == Approx(-0.1724891134).margin(1e-10));
  REQUIRE(Li2(-0.17) == Approx(-0.1632737713).margin(1e-10));
  REQUIRE(Li2(-0.16) == Approx(-0.1540179282).margin(1e-10));
  REQUIRE(Li2(-0.15) == Approx(-0.1447211118).margin(1e-10));
  REQUIRE(Li2(-0.14) == Approx(-0.1353828405).margin(1e-10));
  REQUIRE(Li2(-0.13) == Approx(-0.1260026232).margin(1e-10));
  REQUIRE(Li2(-0.12) == Approx(-0.1165799591).margin(1e-10));
  REQUIRE(Li2(-0.11) == Approx(-0.1071143369).margin(1e-10));
  REQUIRE(Li2(-0.1) == Approx(-0.09760523523).margin(1e-10));
  REQUIRE(Li2(-0.09) == Approx(-0.08805212172).margin(1e-10));
  REQUIRE(Li2(-0.08) == Approx(-0.07845445308).margin(1e-10));
  REQUIRE(Li2(-0.07) == Approx(-0.06881167461).margin(1e-10));
  REQUIRE(Li2(-0.06) == Approx(-0.05912321986).margin(1e-10));
  REQUIRE(Li2(-0.05) == Approx(-0.04938851035).margin(1e-10));
  REQUIRE(Li2(-0.04) == Approx(-0.0396069551).margin(1e-10));
  REQUIRE(Li2(-0.03) == Approx(-0.02977795033).margin(1e-10));
  REQUIRE(Li2(-0.02) == Approx(-0.01990087902).margin(1e-10));
  REQUIRE(Li2(-0.01) == Approx(-0.00997511049).margin(1e-10));
  REQUIRE(Li2(0.) == Approx(0.).margin(1e-10));
  REQUIRE(Li2(0.01) == Approx(0.01002511174).margin(1e-10));
  REQUIRE(Li2(0.02) == Approx(0.02010089902).margin(1e-10));
  REQUIRE(Li2(0.03) == Approx(0.03022805162).margin(1e-10));
  REQUIRE(Li2(0.04) == Approx(0.04040727532).margin(1e-10));
  REQUIRE(Li2(0.05) == Approx(0.05063929246).margin(1e-10));
  REQUIRE(Li2(0.06) == Approx(0.06092484246).margin(1e-10));
  REQUIRE(Li2(0.07) == Approx(0.07126468241).margin(1e-10));
  REQUIRE(Li2(0.08) == Approx(0.0816595877).margin(1e-10));
  REQUIRE(Li2(0.09) == Approx(0.09211035263).margin(1e-10));
  REQUIRE(Li2(0.1) == Approx(0.1026177911).margin(1e-10));
  REQUIRE(Li2(0.11) == Approx(0.1131827373).margin(1e-10));
  REQUIRE(Li2(0.12) == Approx(0.1238060463).margin(1e-10));
  REQUIRE(Li2(0.13) == Approx(0.1344885952).margin(1e-10));
  REQUIRE(Li2(0.14) == Approx(0.1452312834).margin(1e-10));
  REQUIRE(Li2(0.15) == Approx(0.1560350339).margin(1e-10));
  REQUIRE(Li2(0.16) == Approx(0.1669007939).margin(1e-10));
  REQUIRE(Li2(0.17) == Approx(0.1778295358).margin(1e-10));
  REQUIRE(Li2(0.18) == Approx(0.1888222581).margin(1e-10));
  REQUIRE(Li2(0.19) == Approx(0.1998799866).margin(1e-10));
  REQUIRE(Li2(0.2) == Approx(0.2110037754).margin(1e-10));
  REQUIRE(Li2(0.21) == Approx(0.2221947079).margin(1e-10));
  REQUIRE(Li2(0.22) == Approx(0.233453898).margin(1e-10));
  REQUIRE(Li2(0.23) == Approx(0.2447824916).margin(1e-10));
  REQUIRE(Li2(0.24) == Approx(0.2561816675).margin(1e-10));
  REQUIRE(Li2(0.25) == Approx(0.2676526391).margin(1e-10));
  REQUIRE(Li2(0.26) == Approx(0.2791966559).margin(1e-10));
  REQUIRE(Li2(0.27) == Approx(0.2908150047).margin(1e-10));
  REQUIRE(Li2(0.28) == Approx(0.3025090116).margin(1e-10));
  REQUIRE(Li2(0.29) == Approx(0.3142800435).margin(1e-10));
  REQUIRE(Li2(0.3) == Approx(0.3261295101).margin(1e-10));
  REQUIRE(Li2(0.31) == Approx(0.3380588655).margin(1e-10));
  REQUIRE(Li2(0.32) == Approx(0.3500696107).margin(1e-10));
  REQUIRE(Li2(0.33) == Approx(0.3621632955).margin(1e-10));
  REQUIRE(Li2(0.34) == Approx(0.3743415208).margin(1e-10));
  REQUIRE(Li2(0.35) == Approx(0.3866059412).margin(1e-10));
  REQUIRE(Li2(0.36) == Approx(0.3989582673).margin(1e-10));
  REQUIRE(Li2(0.37) == Approx(0.4114002691).margin(1e-10));
  REQUIRE(Li2(0.38) == Approx(0.4239337783).margin(1e-10));
  REQUIRE(Li2(0.39) == Approx(0.4365606916).margin(1e-10));
  REQUIRE(Li2(0.4) == Approx(0.4492829745).margin(1e-10));
  REQUIRE(Li2(0.41) == Approx(0.4621026643).margin(1e-10));
  REQUIRE(Li2(0.42) == Approx(0.4750218745).margin(1e-10));
  REQUIRE(Li2(0.43) == Approx(0.4880427986).margin(1e-10));
  REQUIRE(Li2(0.44) == Approx(0.5011677144).margin(1e-10));
  REQUIRE(Li2(0.45) == Approx(0.5143989892).margin(1e-10));
  REQUIRE(Li2(0.46) == Approx(0.5277390845).margin(1e-10));
  REQUIRE(Li2(0.47) == Approx(0.5411905619).margin(1e-10));
  REQUIRE(Li2(0.48) == Approx(0.5547560886).margin(1e-10));
  REQUIRE(Li2(0.49) == Approx(0.5684384439).margin(1e-10));
  REQUIRE(Li2(0.5) == Approx(0.5822405265).margin(1e-10));
  REQUIRE(Li2(0.51) == Approx(0.5961653614).margin(1e-10));
  REQUIRE(Li2(0.52) == Approx(0.6102161084).margin(1e-10));
  REQUIRE(Li2(0.53) == Approx(0.624396071).margin(1e-10));
  REQUIRE(Li2(0.54) == Approx(0.6387087054).margin(1e-10));
  REQUIRE(Li2(0.55) == Approx(0.6531576315).margin(1e-10));
  REQUIRE(Li2(0.56) == Approx(0.6677466442).margin(1e-10));
  REQUIRE(Li2(0.57) == Approx(0.6824797254).margin(1e-10));
  REQUIRE(Li2(0.58) == Approx(0.6973610584).margin(1e-10));
  REQUIRE(Li2(0.59) == Approx(0.712395042).margin(1e-10));
  REQUIRE(Li2(0.6) == Approx(0.7275863077).margin(1e-10));
  REQUIRE(Li2(0.61) == Approx(0.7429397374).margin(1e-10));
  REQUIRE(Li2(0.62) == Approx(0.7584604836).margin(1e-10));
  REQUIRE(Li2(0.63) == Approx(0.7741539916).margin(1e-10));
  REQUIRE(Li2(0.64) == Approx(0.7900260243).margin(1e-10));
  REQUIRE(Li2(0.65) == Approx(0.8060826895).margin(1e-10));
  REQUIRE(Li2(0.66) == Approx(0.8223304706).margin(1e-10));
  REQUIRE(Li2(0.67) == Approx(0.8387762613).margin(1e-10));
  REQUIRE(Li2(0.68) == Approx(0.8554274037).margin(1e-10));
  REQUIRE(Li2(0.69) == Approx(0.8722917326).margin(1e-10));
  REQUIRE(Li2(0.7) == Approx(0.8893776243).margin(1e-10));
  REQUIRE(Li2(0.71) == Approx(0.9066940527).margin(1e-10));
  REQUIRE(Li2(0.72) == Approx(0.9242506536).margin(1e-10));
  REQUIRE(Li2(0.73) == Approx(0.9420577978).margin(1e-10));
  REQUIRE(Li2(0.74) == Approx(0.9601266752).margin(1e-10));
  REQUIRE(Li2(0.75) == Approx(0.9784693929).margin(1e-10));
  REQUIRE(Li2(0.76) == Approx(0.9970990883).margin(1e-10));
  REQUIRE(Li2(0.77) == Approx(1.016030062).margin(1e-10));
  REQUIRE(Li2(0.78) == Approx(1.035277934).margin(1e-10));
  REQUIRE(Li2(0.79) == Approx(1.05485983).margin(1e-10));
  REQUIRE(Li2(0.8) == Approx(1.0747946).margin(1e-10));
  REQUIRE(Li2(0.81) == Approx(1.095103088).margin(1e-10));
  REQUIRE(Li2(0.82) == Approx(1.115808451).margin(1e-10));
  REQUIRE(Li2(0.83) == Approx(1.13693656).margin(1e-10));
  REQUIRE(Li2(0.84) == Approx(1.158516488).margin(1e-10));
  REQUIRE(Li2(0.85) == Approx(1.180581124).margin(1e-10));
  REQUIRE(Li2(0.86) == Approx(1.203167961).margin(1e-10));
  REQUIRE(Li2(0.87) == Approx(1.226320101).margin(1e-10));
  REQUIRE(Li2(0.88) == Approx(1.250087584).margin(1e-10));
  REQUIRE(Li2(0.89) == Approx(1.27452916).margin(1e-10));
  REQUIRE(Li2(0.9) == Approx(1.299714723).margin(1e-10));
  REQUIRE(Li2(0.91) == Approx(1.325728728).margin(1e-10));
  REQUIRE(Li2(0.92) == Approx(1.352675161).margin(1e-10));
  REQUIRE(Li2(0.93) == Approx(1.380685041).margin(1e-10));
  REQUIRE(Li2(0.94) == Approx(1.4099283).margin(1e-10));
  REQUIRE(Li2(0.95) == Approx(1.440633797).margin(1e-10));
  REQUIRE(Li2(0.96) == Approx(1.47312586).margin(1e-10));
  REQUIRE(Li2(0.97) == Approx(1.507899041).margin(1e-10));
  REQUIRE(Li2(0.98) == Approx(1.545799712).margin(1e-10));
  REQUIRE(Li2(0.99) == Approx(1.588625448).margin(1e-10));
  REQUIRE(Li2(1.) == Approx(1.644934067).margin(1e-10));
}