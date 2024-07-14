import subprocess
import sys
import os
import shutil
from enum import Enum
import fileinput
from argparse import ArgumentParser, ArgumentTypeError
import platform


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


def get_compiler():
    compiler = ""

    if sys.platform != "linux" and sys.platform != "darwin":
        return compiler

    if sys.platform == "linux":
        compiler += "-gcc-"

    if sys.platform == "darwin":
        compiler += "-clang-"
    compiler += get_compiler_version()

    return compiler


def get_profile(os: str, arch: str, build_type: BuildMode):
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

    profile += get_compiler()

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


def get_compiler_version():
    if sys.platform == "linux":
        version_response = subprocess.check_output(
            "gcc --version".split(), encoding="UTF-8"
        ).partition("\n")[0]
        semver_string = version_response[version_response.rfind(" ") + 1 :]
        return semver_string.partition(".")[0]
    if sys.platform == "darwin":
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


def conan_install(profile, additional_options=[], build_missing=False):
    config_settings = [
        "tools.cmake.cmake_layout:build_folder_vars=['settings.os','settings.arch','settings.build_type']"
    ]

    build_profile = get_profile(sys.platform, get_arch(), BuildMode.release)

    cmd = f"conan install . -pr:h BSMPT/{profile} -pr:b BSMPT/{build_profile} ".split()

    for option in additional_options:
        cmd += ["--options", option]

    for conf in config_settings:
        cmd += ["-c", conf]

    if build_missing:
        cmd += ["--build=missing"]

    subprocess.check_call(cmd)


def conan_install_all(
    mode: BuildMode, options=[], build_missing=False, custom_profile=""
):
    if mode == BuildMode.all or mode == BuildMode.release:
        profile = (
            custom_profile
            if custom_profile != ""
            else get_profile(sys.platform, get_arch(), BuildMode.release)
        )
        check_profile(profile)
        conan_install(
            profile,
            options,
            build_missing,
        )
    if mode == BuildMode.all or mode == BuildMode.debug:
        profile = (
            custom_profile
            if custom_profile != ""
            else get_profile(sys.platform, get_arch(), BuildMode.debug)
        )
        check_profile(profile)
        conan_install(
            profile,
            options,
            build_missing,
        )

def create_cmaes():
    cmd = "conan export . --version=0.10.0".split()
    subprocess.check_call(cmd, cwd="tools/conan/cmaes")

def create(build_missing=False):

    config_settings = [
        "tools.cmake.cmake_layout:build_folder_vars=['settings.os','settings.arch','settings.build_type']"
    ]

    profile = get_profile(sys.platform, get_arch(), BuildMode.release)
    cmd = f"conan create . -pr:h BSMPT/{profile} -pr:b BSMPT/{profile}".split()

    for conf in config_settings:
        cmd += ["-c", conf]

    if build_missing:
        cmd += ["--build=missing"]

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


if __name__ == "__main__":

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

    opts = parser.parse_args()

    setup_profiles()

    create_cmaes()

    if opts.create:
        create(build_missing=opts.build_missing,)
    else:

        
        conan_install_all(
            opts.mode,
            opts.options if opts.options is not None else [],
            opts.build_missing,
            opts.profile,
        )
