import subprocess
import sys
import os
import shutil
from enum import Enum
import argparse
from argparse import ArgumentParser


class ArgTypeEnum(Enum):
    @classmethod
    def argtype(cls, s: str) -> Enum:
        try:
            return cls[s]
        except KeyError:
            raise argparse.ArgumentTypeError(f"{s!r} is not a valid {cls.__name__}")

    def __str__(self):
        return self.name


class BuildMode(ArgTypeEnum, Enum):
    all = (0,)
    release = (1,)
    debug = 2


build_profiles = {"win32": "windows-release", "linux": "linux-release", "darwin": "macos-release"}

target_profiles = {
    "win32": {"release": "windows-release", "debug": "windows-debug"},
    "linux": {"release": "linux-release", "debug": "linux-debug"},
    "darwin": {"release": "macos-release", "debug": "macos-debug"}
}


def setup_profiles():
    conan_home = subprocess.check_output("conan config home".split(), encoding="UTF-8")
    profile_dir = os.path.join(str(conan_home.split()[0]), "profiles", "BSMPT")
    print(profile_dir)
    try:
        shutil.copytree("profiles/BSMPT", profile_dir)
    except:
        pass


def conan_install(profile, additional_options=[]):
    config_settings = [
        "tools.cmake.cmake_layout:build_folder_vars=['settings.os','settings.arch','settings.build_type']"
    ]

    build_profile = build_profiles[sys.platform]

    cmd = f"conan install . -pr:h BSMPT/{profile} -pr:b BSMPT/{build_profile} ".split()

    for option in additional_options:
        cmd += ["--options", option]

    for conf in config_settings:
        cmd += ["-c", conf]

    subprocess.check_call(cmd)


def conan_install_all(mode: BuildMode, options=[]):
    profiles = target_profiles[sys.platform]
    if mode == BuildMode.all or mode == BuildMode.release:
        conan_install(profiles["release"],options)
    if mode == BuildMode.all or mode == BuildMode.debug:
        conan_install(profiles["debug"],options)


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
        "--mode", default=BuildMode.release, type=BuildMode.argtype, choices=BuildMode
    )
    parser.add_argument("--options", nargs='+')

    opts = parser.parse_args()
    setup_profiles()
    options=opts.options
    if options is None:
        options = []
    conan_install_all(opts.mode, opts.options)
