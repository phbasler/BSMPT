from conan import ConanFile
from conan.tools.cmake import cmake_layout

class BSMPT(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps"

    options = { "enable_tests": [True, False]}
    default_options = {"enable_tests": True}
    
    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires("boost/1.84.0")
        self.requires("gsl/2.7.1")
        self.requires("nlohmann_json/3.11.3")
        self.requires("nlopt/2.7.1")

    def build_requirements(self):
        self.tool_requires("cmake/3.29.0")

    def layout(self):
        cmake_layout(self)