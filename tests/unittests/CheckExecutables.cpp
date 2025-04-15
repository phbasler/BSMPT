#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <vector>
#include <iostream>

namespace fs = std::filesystem;

TEST_CASE("Check for existence of BSMPT executable", "[executable]")
{

  std::vector<fs::path> matching_files;

  // Iterate through directories inside 'build'
  fs::path base_directory = BASE_PATH;
  base_directory /= "build";
  std::cout << "bp = " << base_directory << std::endl;
  for (const auto &dir_entry : fs::directory_iterator(base_directory))
  {
    if (fs::is_directory(dir_entry.path()))
    {
      fs::path file_path = dir_entry.path() / "bin" / "BSMPT";
      if (fs::exists(file_path) && fs::is_regular_file(file_path))
      {
        matching_files.push_back(file_path);
      }
    }
  }

  REQUIRE_FALSE(matching_files.empty());
}