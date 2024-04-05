import subprocess


def conan_install():
    config_settings = ["tools.cmake.cmake_layout:build_folder_vars=['settings.os','settings.arch','settings.build_type']"]
    cmd = f"conan install . ".split()
    for conf in config_settings:
        cmd += ["-c", conf]

    subprocess.check_call(cmd)

if __name__ == "__main__":
    conan_install()