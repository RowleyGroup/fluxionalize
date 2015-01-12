#!/usr/bin/python
import os, sys

#sets up system
def set_up(config_file):
        set_path = open('set.sh', 'w')
        set_path.write("#!/bin/bash\n")
        if config_file:
                if not os.path.exists(config_file):
                        print("Invalid config file!")
                        return 0

                config = open(config_file, 'r')
                for line in config:
                        prog = line.split()
                        if is_exe(prog[1]):
                                set_path.write('export PATH='+prog[1]+':$PATH\n')
                                print("is executable!")
        else:
                if not path_finder("obabel"):
                        print("Open Babel not found!")
                if not path_finder("antechamber"):
                        print("Antechamber not found!")
                vmd_path=path_finder("vmd")
                if not vmd_path:
                        print("VMD not found!")
                elif vmd_path:
                        if not path_finder("catdcd"):
                                set_path.write('export PATH='+vmd_path[:-9]+'/vmd-1.9.1/plugins/LINUXAMD64/bin/catdcd4.0:$PATH')
                if not path_finder("namd2"):
                        print("NAMD not found!")
        set_path.close()

#checks if executables in path, adds them if they're not there
def path_finder(program):
        for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                        print ("Executable in PATH as %s!")%exe_file
                        return exe_file
        return 0

#checks if given path is an existing executable file
def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

infile = sys.argv[1]

set_up(infile)
