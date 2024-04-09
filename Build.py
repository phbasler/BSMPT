import Setup
import subprocess
import sys

presets={
    "linux": "conan-linux-x86_64-release",
    "win32": "conan-windows-x86_64-release",
    "darwin": "conan-macos-x86_64-release"
}

def build(preset):
    cmd=f"cmake --preset {preset}".split()
    subprocess.check_call(cmd)

    cmd=f"cmake --build --preset {preset}".split()
    subprocess.check_call(cmd)


def main():
    Setup.conan_install_all(Setup.BuildMode.release)
    build(presets[sys.platform])

if __name__ == "__main__":
    main()
    