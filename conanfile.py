from conan import ConanFile, tools
from conan.tools.system.package_manager import Apt
from conan.errors import ConanInvalidConfiguration
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conan.tools.files import load
from conan.tools.scm import Git
import os, re
from conan.tools.env import Environment
from conan.errors import ConanException

required_conan_version = ">=2.0.0 <3"


class BSMPT(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps"

    name = "bsmpt"

    exports_sources = (
        "CMakeLists.txt",
        "src/*",
        "include/*",
        "tools/*",
        "tests/*",
        "standalone/*",
    )

    options = {
        "fPIC": [True, False],
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
        "UseVectorization": [True, False],  # use vectorization for the build,
        "BuildExecutables": [True, False],
    }
    default_options = {
        "fPIC": True,
        "EnableTests": True,
        "UseLibCMAES": True,
        "UseNLopt": True,
        "MakeAdditionalTesting": False,
        "CompileBaryo": False,
        "EnableCoverage": False,
        "UseVectorization": False,  # This causes double free with cmaes. We need to modify the package to include vectorization as well as an option
        "BuildExecutables": True,
    }

    def requirements(self):
        self.requires("eigen/3.4.0", transitive_headers=True, transitive_libs=True)
        self.requires("gsl/2.7.1", transitive_headers=True, transitive_libs=True)
        self.requires("nlohmann_json/3.11.3", transitive_headers=True)

        if self.options.CompileBaryo:
            self.requires("boost/1.84.0", transitive_headers=True, transitive_libs=True)

        if self.options.UseNLopt:
            self.requires("nlopt/2.9.1", transitive_headers=True, transitive_libs=True)

        if self.options.UseLibCMAES:
            self.requires(
                "cmaes/0.10.0@bsmpt/local",
                transitive_headers=True,
                transitive_libs=True,
            )

    def build_requirements(self):
        self.tool_requires("cmake/3.29.0")

        if self.options.EnableTests:
            self.test_requires("catch2/3.5.3")
            self.test_requires("benchmark/1.6.2")

    def system_requirements(self):
        if self.options.EnableCoverage:
            apt = Apt(self)
            apt.install(["lcov"], update=True, check=True)

    def set_version(self):
        content = load(self, os.path.join(self.recipe_folder, "CMakeLists.txt"))
        value = re.search(r"set\(BSMPT_VERSION (.*)\)", content)
        extracted_version = value.group(1).strip()

        is_git_tag = False
        git = Git(self, folder=self.recipe_folder)
        try:
            git.run("describe --exact-match --tags")
            is_git_tag = True
        except Exception:
            is_git_tag = False

        if is_git_tag:
            self.version = extracted_version
        else:
            # if not tag -> pre-release version
            try:
                commit_hash = git.get_commit()[:8]
                self.version = f"{extracted_version}.{commit_hash}"
            except ConanException:
                # In this case (no git tag but also no git available) the source code was downloaded in a different way.
                # We don't know if it is a changed code or the zip from the download, so we stick to the cmake defined version
                self.version = extracted_version

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
        tc.variables["BSMPTBuildExecutables"] = self.options.BuildExecutables

        tc.generate()

    def validate(self):
        if self.options.UseVectorization and self.options.EnableCoverage:
            raise ConanInvalidConfiguration(
                "Vectorization and coverage are not supported simultaneously."
            )

        if self.settings.os != "Linux" and self.options.EnableCoverage:
            raise ConanInvalidConfiguration("We depend on lcov for coverage.")

        if (
            self.settings.os == "Linux"
            and self.options.UseLibCMAES
            and self.options.UseVectorization
        ):
            raise ConanInvalidConfiguration(
                "This causes a double free error. CMAES needs to be modified to be build with Vectorization as well"
            )

        tools.build.check_min_cppstd(self, "17")

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.parallel = True
        cmake.build()

        environment = Environment()
        environment.define("CTEST_OUTPUT_ON_FAILURE", "1")
        environment.define("CTEST_PARALLEL_LEVEL", str(os.cpu_count()))
        envvars = environment.vars(self)

        if self.options.get_safe("EnableTests"):
            with envvars.apply():
                cmake.test()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_id(self):
        del self.info.options.EnableTests

    def package_info(self):
        self.cpp_info.components["ASCIIPlotter"].libs = ["ASCIIPlotter"]
        self.cpp_info.components["ASCIIPlotter"].requires = [
            "nlohmann_json::nlohmann_json",
        ]
        self.cpp_info.components["ASCIIPlotter"].set_property(
            "cmake_target_name", "BSMPT::ASCIIPlotter"
        )

        if self.options.CompileBaryo:
            self.cpp_info.components["ASCIIPlotter"].requires.append(
                "boost::boost",
            )

        self.cpp_info.components["Spline"].libs = ["Spline"]
        self.cpp_info.components["Spline"].requires = [
            "nlohmann_json::nlohmann_json",
        ]
        self.cpp_info.components["Spline"].set_property(
            "cmake_target_name", "BSMPT::Spline"
        )

        if self.options.CompileBaryo:
            self.cpp_info.components["Spline"].requires.append(
                "boost::boost",
            )

        self.cpp_info.components["Utility"].libs = ["Utility"]
        self.cpp_info.components["Utility"].requires = [
            "gsl::gsl",
            "nlohmann_json::nlohmann_json",
            "ASCIIPlotter",
            "Spline",
        ]
        self.cpp_info.components["Utility"].set_property(
            "cmake_target_name", "BSMPT::Utility"
        )

        if self.options.CompileBaryo:
            self.cpp_info.components["Utility"].requires.append(
                "boost::boost",
            )

        self.cpp_info.components["BounceSolution"].libs = ["BounceSolution"]
        self.cpp_info.components["BounceSolution"].requires = [
            "eigen::eigen",
            "gsl::gsl",
            "Minimizer",
            "Utility",
            "MinimumTracer",
        ]
        self.cpp_info.components["BounceSolution"].set_property(
            "cmake_target_name", "BSMPT::BounceSolution"
        )

        self.cpp_info.components["GW"].libs = ["GW"]
        self.cpp_info.components["GW"].requires = [
            "eigen::eigen",
            "gsl::gsl",
            "Minimizer",
            "Utility",
            "BounceSolution",
        ]
        self.cpp_info.components["GW"].set_property("cmake_target_name", "BSMPT::GW")

        self.cpp_info.components["Minimizer"].libs = ["Minimizer"]
        self.cpp_info.components["Minimizer"].requires = [
            "eigen::eigen",
            "gsl::gsl",
            # "Threads::Threads",
            "Utility",
            # "Models",
        ]
        self.cpp_info.components["Minimizer"].set_property(
            "cmake_target_name", "BSMPT::Minimizer"
        )

        if self.options.UseNLopt:

            self.cpp_info.components["Minimizer"].requires.append("nlopt::nlopt")

        if self.options.UseLibCMAES:
            self.cpp_info.components["Minimizer"].requires.append("cmaes::cmaes")

        self.cpp_info.components["MinimumTracer"].libs = ["MinimumTracer"]
        self.cpp_info.components["MinimumTracer"].requires = [
            "eigen::eigen",
            "gsl::gsl",
            "Minimizer",
            "Utility",
        ]
        self.cpp_info.components["MinimumTracer"].set_property(
            "cmake_target_name", "BSMPT::MinimumTracer"
        )

        self.cpp_info.components["Models"].libs = ["Models"]
        self.cpp_info.components["Models"].requires = [
            "gsl::gsl",
            "eigen::eigen",
            "Minimizer",
            "ThermalFunctions",
            "Utility",
        ]
        self.cpp_info.components["Models"].set_property(
            "cmake_target_name", "BSMPT::Models"
        )

        self.cpp_info.components["ThermalFunctions"].libs = ["ThermalFunctions"]
        self.cpp_info.components["ThermalFunctions"].requires = ["gsl::gsl", "Utility"]
        self.cpp_info.components["ThermalFunctions"].set_property(
            "cmake_target_name", "BSMPT::ThermalFunctions"
        )

        self.cpp_info.components["TransitionTracer"].libs = ["TransitionTracer"]
        self.cpp_info.components["TransitionTracer"].requires = [
            "BounceSolution",
            "MinimumTracer",
            "GW",
        ]
        self.cpp_info.components["TransitionTracer"].set_property(
            "cmake_target_name", "BSMPT::TransitionTracer"
        )
