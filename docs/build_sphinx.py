import os
import sys
from sphinx.cmd.build import main as sphinx_main

def main():
    if len(sys.argv) < 2:
        print("Usage: python build.py /path/to/xml [sphinx options]") 
        sys.exit(1)
    # Default values for Sphinx options

    # Get the directory where the script is located 
    script_dir = os.path.dirname(os.path.abspath(__file__))

    defaults = {
        'builder': 'html',  # Default builder
        'config':  script_dir,  # Default configuration directory
        'source': script_dir,  # Default source directory
        'build': os.path.join(script_dir,"_build")  # Default build directory
    }

    # Convert default XML directory to an absolute path
    xml_dir = os.path.abspath(sys.argv[1])
    os.environ['XML_DIR'] = xml_dir

    # Collect command-line arguments
    cli_args = sys.argv[2:]

    # Update defaults with command-line arguments if provided
    if '-b' in cli_args:
        defaults['builder'] = cli_args[cli_args.index('-b') + 1]
    if '-c' in cli_args:
        defaults['config'] = cli_args[cli_args.index('-c') + 1]
    if len(cli_args) > 0:
        defaults['source'] = cli_args[0]
    if len(cli_args) > 1:
        defaults['build'] = cli_args[1]

     # Convert source and build paths to absolute paths 
    defaults['config'] = os.path.abspath(defaults['config']) 
    defaults['source'] = os.path.abspath(defaults['source']) 
    defaults['build'] = os.path.abspath(defaults['build'])

    # Construct the Sphinx command
    sphinx_cmd = [
        # 'sphinx-build', 
        '-b', defaults['builder'], 
        '-c', defaults['config'],
        defaults['source'],
        defaults['build']
    ]

    print(sphinx_cmd)

    # Run Sphinx build with the constructed command
    sphinx_main(sphinx_cmd)

if __name__ == "__main__":
    main()
