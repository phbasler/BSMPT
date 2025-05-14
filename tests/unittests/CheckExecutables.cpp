#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <iostream>
#include <vector>

namespace fs = std::filesystem;

TEST_CASE("Check for existence of BSMPT executable", "[executable]")
{

  std::vector<fs::path> matching_files;

  std::cout << "Start test" << std::endl;

  // Iterate through directories inside 'build'
  fs::path base_directory = BASE_PATH;
  base_directory /= "build";
  std::cout << "bp = " << base_directory << std::endl;
  for (const auto &dir_entry : fs::directory_iterator(base_directory))
  {
    std::cout << "Checking entry: " << dir_entry.path() << std::endl;
    if (fs::is_directory(dir_entry.path()))
    {
      fs::path file_path = dir_entry.path() / "bin";
      // Append ".exe" only on Windows

      for (const auto &bin_dir : fs::directory_iterator(file_path))
      {
        fs::path low_file_path = bin_dir.path() / "BSMPT";
#ifdef _WIN32
        low_file_path.replace_extension(".exe");
#endif
        if (fs::exists(low_file_path) and fs::is_regular_file(low_file_path))
        {
          std::cout << low_file_path << " is a regular file" << std::endl;
          matching_files.push_back(low_file_path);
        }
      }

      file_path /= "BSMPT";
#ifdef _WIN32
      file_path.replace_extension(".exe");
#endif
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