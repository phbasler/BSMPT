// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this
                          // in one cpp file
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include "catch.hpp"
#include <BSMPT/utility/Logger.h>
struct MyListener : Catch::TestEventListenerBase
{

  using TestEventListenerBase::TestEventListenerBase; // inherit constructor

  void testCaseStarting(Catch::TestCaseInfo const &testInfo) override
  {
    (void)testInfo;
    BSMPT::Logger::SetLevel(BSMPT::LoggingLevel::ProgDetailed, true);
  }

  //  void testCaseEnded(Catch::TestCaseStats const &testCaseStats) override
  //  {
  //    // Tear-down after a test case is run
  //  }
};
CATCH_REGISTER_LISTENER(MyListener)
