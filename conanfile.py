from conan import ConanFile, tools
from conan.tools.system.package_manager import Apt
from conan.errors import ConanInvalidConfiguration
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conan.tools.files import load, update_conandata
from conan.tools.scm import Git

required_conan_version = ">=2.0.0 <3"


class BSMPT(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps"

    name = "bsmpt"
    version = "3.0.2"

    exports_sources = "CMakeLists.txt", "src/*", "include/*", "tools/*", "tests/*" , "standalone/*"

    options = {
        "EnableTests": [True, False],  # enables the unit tests
        "UseLibCMAES": [
            True,
            False,
        ],  # Use CMAES. Fetches it through cmake_fetch if not installed
        "UseNLopt": [True, False],  # Use NLopt for minimization
        "MakeAdditionalTesting": [True, False],  # build additional test executables
        "CompileBaryo": [
            True,
            False,
        ],  # compile the electroweak baryogenesis for the C2HDM
        "EnableCoverage": [True, False],  # enable code coverage
        "UseVectorization": [True, False],  # use vectorization for the build
    }
    default_options = {
        "EnableTests": True,
        "UseLibCMAES": True,
        "UseNLopt": True,
        "MakeAdditionalTesting": False,
        "CompileBaryo": False,
        "EnableCoverage": False,
        "UseVectorization": True,
    }

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires("gsl/2.7.1")
        self.requires("nlohmann_json/3.11.3")

        if self.options.CompileBaryo:
            self.requires("boost/1.84.0")

        if self.options.UseNLopt:
            self.requires("nlopt/2.7.1")

    def build_requirements(self):
        self.tool_requires("cmake/3.29.0")

        if self.options.EnableTests:
            self.test_requires("catch2/3.5.3")
            self.test_requires("benchmark/1.6.2")

    def system_requirements(self):
        if self.options.EnableCoverage:
            apt = Apt(self)
            apt.install(["lcov"], update=True, check=True)

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)

        tc.variables["BUILD_TESTING"] = self.options.EnableTests
        tc.variables["UseLibCMAES"] = self.options.UseLibCMAES
        tc.variables["UseNLopt"] = self.options.UseNLopt
        tc.variables["MakeAdditionalTesting"] = self.options.MakeAdditionalTesting
        tc.variables["BSMPTCompileBaryo"] = self.options.CompileBaryo
        tc.variables["EnableCoverage"] = self.options.EnableCoverage
        tc.variables["BSMPTUseVectorization"] = self.options.UseVectorization

        tc.generate()

    def validate(self):
        if self.options.UseVectorization and self.options.EnableCoverage:
            raise ConanInvalidConfiguration(
                "Vectorization and coverage are not supported simultaneously."
            )

        if self.settings.os != "Linux" and self.options.EnableCoverage:
            raise ConanInvalidConfiguration("We depend on lcov for coverage.")

        tools.build.check_min_cppstd(self, "17")


    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC



    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()


