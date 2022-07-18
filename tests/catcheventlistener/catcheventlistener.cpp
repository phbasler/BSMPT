// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/utility/Logger.h>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
struct MyListener : Catch::EventListenerBase
{

  using EventListenerBase::EventListenerBase; // inherit constructor

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
