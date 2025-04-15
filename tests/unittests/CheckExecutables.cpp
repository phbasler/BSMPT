#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <iostream>
#include <vector>

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
    std::cout << "Checking entry: " << dir_entry.path() << std::endl;
    if (fs::is_directory(dir_entry.path()))
    {
      fs::path file_path = dir_entry.path() / "bin" / "BSMPT";
      std::cout << "Checking file" << file_path << std::endl;
      if (fs::exists(file_path))
      {
        std::cout << "File exists" << std::endl;
        if (fs::is_regular_file(file_path))
        {
          std::cout << file_path << " is a regular file" << std::endl;
          matching_files.push_back(file_path);
        }
      }
    }
  }

  REQUIRE_FALSE(matching_files.empty());
}