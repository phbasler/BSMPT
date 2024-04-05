from conan import ConanFile

class BSMPT(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps"
    
    def requirements(self):
        self.requires("eigen/3.4.0")