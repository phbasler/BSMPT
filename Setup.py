import subprocess
import sys
import os
import shutil
from enum import Enum
import fileinput
from argparse import ArgumentParser, ArgumentTypeError
import platform
from pathlib import Path


class ArgTypeEnum(Enum):
    @classmethod
    def argtype(cls, s: str) -> Enum:
        try:
            return cls[s]
        except KeyError:
            raise ArgumentTypeError(f"{s!r} is not a valid {cls.__name__}")

    def __str__(self):
        return self.name


class BuildMode(ArgTypeEnum, Enum):
    all = (0,)
    release = (1,)
    debug = 2


class Compiler(ArgTypeEnum, Enum):
    gcc = (0,)
    clang = 1


def get_compiler(compiler: Compiler):
    compilerString = ""

    if sys.platform != "linux" and sys.platform != "darwin":
        return compilerString

    if sys.platform == "linux":
        if compiler == Compiler.gcc or compiler is None:
            compilerString += "-gcc-"
        else:
            compilerString += "-clang-"

    if sys.platform == "darwin":
        compilerString += "-clang-"
    compilerString += get_compiler_version(compiler)

    return compilerString


def get_profile(os: str, arch: str, build_type: BuildMode, compiler: Compiler):
    profile = ""
    if os == "win32":
        profile += "windows"
    elif os == "linux":
        profile += "linux"
    elif os == "darwin":
        profile += "macos"

    profile += "-"

    if build_type == BuildMode.release:
        profile += "release"
    elif build_type == BuildMode.debug:
        profile += "debug"

    profile += "-"
    profile += arch

    profile += get_compiler(compiler)

    return profile


def set_setting(file, setting, value):
    for line in fileinput.input([file], inplace=True):
        if line.strip().startswith(setting):
            line = setting + "=" + value + "\n"
        sys.stdout.write(line)


def check_profile(profile):
    path = os.path.join("profiles", "BSMPT", profile)
    if not os.path.isfile(path):
        conan_home = subprocess.check_output(
            "conan config home".split(), encoding="UTF-8"
        ).split("\n")[0]
        print(
            f"Profile does not exist in BSMPT/profiles.\nUsing profile {profile} created from the default profile. Change it accordingly."
        )
        if not os.path.isfile(conan_home + "/profiles/default"):
            cmd = "conan profile detect".split()
            subprocess.check_output(cmd)
        if sys.platform != "win32":
            cmd = (
                "cp " + conan_home + "/profiles/default profiles/BSMPT/" + str(profile)
            )
            subprocess.check_call(cmd, shell=True)
            set_setting(path, "compiler.cppstd", "gnu17")

        else:
            cmd = (
                "copy "
                + conan_home
                + "\\profiles\\default profiles\\BSMPT\\"
                + str(profile)
            )
            subprocess.check_call(cmd, shell=True)
            set_setting(path, "compiler.cppstd", "17")

        setup_profiles()
        check_profile(profile)


def get_compiler_version(compiler: Compiler):
    if sys.platform == "linux" and compiler != Compiler.clang:
        version_response = subprocess.check_output(
            "gcc --version".split(), encoding="UTF-8"
        ).partition("\n")[0]
        semver_string = version_response[version_response.rfind(" ") + 1 :]
        return semver_string.partition(".")[0]
    elif sys.platform == "darwin" or compiler == Compiler.clang:
        version_response = subprocess.check_output(
            "clang --version".split(), encoding="UTF-8"
        ).partition("\n")[0]
        semver_string = version_response[version_response.rfind("version") + 8 :]
        return semver_string.partition(".")[0]
    return ""


def get_arch():
    arch = "x86_64"
    if platform.machine() == "aarch64" or platform.machine() == "arm64":
        arch = "armv8"
    return arch


def setup_profiles():
    conan_home = subprocess.check_output("conan config home".split(), encoding="UTF-8")
    profile_dir = os.path.join(str(conan_home.split()[0]), "profiles", "BSMPT")
    print(profile_dir)
    if os.path.exists(profile_dir):
        shutil.rmtree(profile_dir)
    shutil.copytree("profiles/BSMPT", profile_dir)

def setup_cmaes():
    file_directory = Path(__file__).parent.absolute()
    cmaes_dir = os.path.join(file_directory, "tools", "conan", "cmaes","all")

    # Define the recipe name
    recipe = "cmaes/0.10.0@bsmpt/local"

    try:
        # Run the conan search command and capture the output
        
        result = subprocess.check_output(f"conan list {recipe} -c".split(), stderr=subprocess.STDOUT, text=True)
        
        
        # Check if the output indicates the recipe is not found
        if "ERROR: Recipe" in result:
            print(f"Recipe '{recipe}' not found. Exporting...")
            subprocess.check_output("conan export conanfile.py --version 0.10.0 --user bsmpt --channel local".split(), cwd=cmaes_dir)
            print(f"Recipe '{recipe}' successfully exported.")
        else:
            print(f"Recipe '{recipe}' already exists in the local cache.")
    except subprocess.CalledProcessError as e:
        # Handle errors from the subprocess
        error_output = e.output.decode("utf-8")  # Decode the error output for debugging
        if "ERROR: Recipe" in error_output:
            print(f"Recipe '{recipe}' not found in local cache. Exporting...")
            subprocess.check_output("conan export conanfile.py --version 0.10.0 --user bsmpt --channel local".split(), cwd=cmaes_dir)
            print(f"Recipe '{recipe}' successfully exported.")
        else:
            # If another error occurs, print the error output
            print(f"An error occurred: {error_output}")

    

    



def conan_install(
    profile, additional_options=[], build_missing=False, compiler: Compiler = None
):
    config_settings = [
        "tools.cmake.cmake_layout:build_folder_vars=['settings.os','settings.arch','settings.build_type']"
    ]

    build_profile = get_profile(sys.platform, get_arch(), BuildMode.release, compiler)

    cmd = f"conan install . -pr:h BSMPT/{profile} -pr:b BSMPT/{build_profile} ".split()

    for option in additional_options:
        cmd += ["--options", option]

    for conf in config_settings:
        cmd += ["-c", conf]

    if build_missing:
        cmd += ["--build=missing"]

    print(f"Executing command {cmd}")

    subprocess.check_call(cmd)


def conan_install_all(
    mode: BuildMode,
    options=[],
    build_missing=False,
    custom_profile="",
    compiler: Compiler = None,
):
    if mode == BuildMode.all or mode == BuildMode.release:
        profile = (
            custom_profile
            if custom_profile != ""
            else get_profile(sys.platform, get_arch(), BuildMode.release, compiler)
        )
        check_profile(profile)
        conan_install(profile, options, build_missing, compiler)
    if mode == BuildMode.all or mode == BuildMode.debug:
        profile = (
            custom_profile
            if custom_profile != ""
            else get_profile(sys.platform, get_arch(), BuildMode.debug, compiler)
        )
        check_profile(profile)
        conan_install(profile, options, build_missing, compiler)


def create(build_missing=False, compiler: Compiler = None, additional_options=[]):

    config_settings = [
        "tools.cmake.cmake_layout:build_folder_vars=['settings.os','settings.arch','settings.build_type']"
    ]

    profile = get_profile(sys.platform, get_arch(), BuildMode.release, compiler)
    cmd = f"conan create . -pr:h BSMPT/{profile} -pr:b BSMPT/{profile}".split()

    for conf in config_settings:
        cmd += ["-c", conf]

    if build_missing:
        cmd += ["--build=missing"]
    

    for option in additional_options:
        cmd += ["--options", option]

    if "EnableTests=True" not in additional_options:
        cmd += ["--options", "EnableTests=False"]
        
    cmd += ["--options", "BuildExecutables=False"]

    subprocess.check_call(cmd)


class ArgTypeEnum(Enum):
    @classmethod
    def argtype(cls, s: str) -> Enum:
        try:
            return cls[s]
        except KeyError:
            raise argparse.ArgumentTypeError(f"{s!r} is not a valid {cls.__name__}")

    def __str__(self):
        return self.name


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument(
        "--mode",
        "-m",
        default=BuildMode.release,
        type=BuildMode.argtype,
        choices=BuildMode,
        help="Should BSMPT be build in Debug, Release or both?",
    )
    parser.add_argument(
        "--options",
        "-o",
        nargs="+",
        action="extend",
        help="Options to pass through to conan. For the available options please look into the conanfile.",
    )
    parser.add_argument(
        "--build-missing", "-b", action="store_true", help="Build missing dependencies."
    )

    parser.add_argument(
        "-p",
        "--profile",
        type=str,
        help="The name of a custom profile. If you leave this empty we will try to deduce a matching profile.",
        default="",
    )

    parser.add_argument(
        "-c", "--create", action="store_true", help="create the local conan package"
    )

    parser.add_argument(
        "-co",
        "--compiler",
        default=None,
        type=Compiler.argtype,
        choices=Compiler,
        help="Force a certain compiler",
    )

    return parser.parse_args()


if __name__ == "__main__":

    opts = parse_arguments()
    setup_cmaes()
    setup_profiles()

    options = []
    if opts.options is not None:
        options = opts.options

    if opts.create:
        create(
            build_missing=opts.build_missing,
            compiler=opts.compiler,
            additional_options=options,
        )
    else:

        conan_install_all(
            opts.mode,
            options,
            opts.build_missing,
            opts.profile,
            opts.compiler,
        )
