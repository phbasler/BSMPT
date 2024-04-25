import os
import sys
 
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import Build

if __name__ == "__main__":
    print(Build.get_preset())