import subprocess
import sys
import os
import shutil
from enum import Enum
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

    return profile


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


def conan_install_all(mode: BuildMode, options=[], build_missing=False):
    if mode == BuildMode.all or mode == BuildMode.release:
        conan_install( get_profile(sys.platform, get_arch(), BuildMode.release) , options, build_missing)
    if mode == BuildMode.all or mode == BuildMode.debug:
        conan_install(get_profile(sys.platform, get_arch(), BuildMode.debug), options, build_missing)


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

    opts = parser.parse_args()
    setup_profiles()
    conan_install_all(
        opts.mode, opts.options if opts.options is not None else [], opts.build_missing
    )
