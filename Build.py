import Setup
import subprocess
import sys


def get_preset(mode : Setup.BuildMode):
    preset = "conan-"
    os = sys.platform
    if os == "win32":
        preset += "windows"
    elif os == "linux":
        preset += "linux"
    elif os == "darwin":
        preset += "macos"

    preset += "-"
    preset += Setup.get_arch()

    if mode == Setup.BuildMode.release:
        preset += "-release"
    elif mode == Setup.BuildMode.debug:
        preset += "-debug"

    return preset


def build(preset):
    cmd=f"cmake --preset {preset} --fresh".split()
    subprocess.check_call(cmd)

    cmd = f"cmake --build --preset {preset}".split()
    subprocess.check_call(cmd)


def main():
    opts = Setup.parse_arguments()
    Setup.setup_profiles()
    Setup.conan_install_all(Setup.BuildMode.release,
                            opts.options if opts.options is not None else [],
                            build_missing=True,
                            custom_profile=opts.profile
                            )
    build(get_preset(Setup.BuildMode.release))


if __name__ == "__main__":
    main()
