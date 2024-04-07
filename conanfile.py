from conan import ConanFile
from conan.tools.cmake import cmake_layout, CMakeToolchain


class BSMPT(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps"

    options = {
        "enable_tests": [True, False],
        "UseLibCMAES": [True, False],
        "UseNLopt": [True, False],
        "MakeAdditionalTesting": [True, False],
        "BSMPTCompileBaryo": [True, False],
        "EnableCoverage": [True, False],
    }
    default_options = {
        "enable_tests": True,
        "UseLibCMAES": True,
        "UseNLopt": True,
        "MakeAdditionalTesting": False,
        "BSMPTCompileBaryo": True,
        "EnableCoverage": True,
    }

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires("boost/1.84.0")
        self.requires("gsl/2.7.1")
        self.requires("nlohmann_json/3.11.3")

        if self.options.UseNLopt:
            self.requires("nlopt/2.7.1")

    def build_requirements(self):
        self.tool_requires("cmake/3.29.0")
        
        if self.options.enable_tests:
            self.test_requires("catch2/3.5.3")
            self.test_requires("benchmark/1.6.1")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)

        tc.variables["BUILD_TESTING"] = self.options.enable_tests
        tc.variables["UseLibCMAES"] = self.options.UseLibCMAES
        tc.variables["UseNLopt"] = self.options.UseNLopt
        tc.variables["MakeAdditionalTesting"] = self.options.MakeAdditionalTesting
        tc.variables["BSMPTCompileBaryo"] = self.options.BSMPTCompileBaryo
        tc.variables["EnableCoverage"] = self.options.EnableCoverage

        tc.generate()
