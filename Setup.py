import subprocess
import sys
import os
import shutil

build_profiles = {
    "win32": "windows-release",
    "linux": "linux-release"
}

target_profiles = {
    "win32": [ "windows-release", "windows-debug"],
    "linux": ["linux-release", "linux-debug"]
}


def setup_profiles():
    conan_home=subprocess.check_output('conan config home'.split(),encoding='UTF-8')
    profile_dir=os.path.join(str(conan_home.split()[0]),"profiles","BSMPT")
    print(profile_dir)
    try:
        shutil.copytree("profiles/BSMPT",profile_dir)
    except:
        pass


def conan_install(profile, additional_options = []):
    config_settings = ["tools.cmake.cmake_layout:build_folder_vars=['settings.os','settings.arch','settings.build_type']"]

    build_profile=build_profiles[sys.platform]

    cmd = f"conan install . -pr:h BSMPT/{profile} -pr:b BSMPT/{build_profile} ".split()

    for option in additional_options:
        cmd += ["--options", option]

    for conf in config_settings:
        cmd += ["-c", conf]

    subprocess.check_call(cmd)

def conan_install_all():
    for profile in target_profiles[sys.platform]:
        conan_install(profile)

if __name__ == "__main__":
    setup_profiles()
    conan_install_all()