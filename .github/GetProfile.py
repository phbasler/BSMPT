import os
import sys
 
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import Build

if __name__ == "__main__":

    env_file = os.getenv('GITHUB_ENV')
    with open(env_file, "a") as myfile:
        myfile.write(f"generated_cmake_preset={Build.get_preset()}")
