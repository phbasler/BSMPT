import Setup
import subprocess
import sys

def get_preset():
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

    preset += "-release"

    return preset

def build(preset):
    cmd=f"cmake --preset {preset}".split()
    subprocess.check_call(cmd)

    cmd=f"cmake --build --preset {preset}".split()
    subprocess.check_call(cmd)


def main():
    Setup.setup_profiles()
    Setup.conan_install_all(Setup.BuildMode.release,build_missing=True)
    build(get_preset())

if __name__ == "__main__":
    main()
    