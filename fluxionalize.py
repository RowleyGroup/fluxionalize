import sys, getopt, os, subprocess

#write test.sh and send it to queue
def genSubScript():
        location = os.getcwd()
        name_spec=False
        num_procs='2'
        num_clusters='1'
        (opts, args)= getopt.getopt(sys.argv[1:], "p:n:l:i:c:")
        file_ext='.pdb'
        
        for opt, arg in opts:
                if opt=='-p':
                        num_procs=arg
                elif opt=='-c':
                        num_clusters=arg
                elif opt=='-n':
                        name=arg
                        name_spec=True
                elif opt=='-l':
                        if os.path.exists(arg)and os.access(arg, os.X_OK):
                                location=arg
                elif opt=='-i':
                        input = arg
                        file_formats=['.mol','.mdl','.sdf','.sd','.pdb','.ent','.ml2','.sy2','.mol2','.xml','.dat','.output','.out','.log','.xyz','.arc','.car','.msi','.pc','.pos','.g09','.g03','.g98','.g94',
                                        '.g92','.gal','.gzmat','.mpo','.xsf','.cdx']

                        (file_name, file_ext)=os.path.splitext(input)

                        if not any (file_ext in type for type in file_formats):
                                print("File type %s may not be supported.") %file_ext
                        if file_ext=="":
                                file_ext='.SMILES'

        if name_spec==False:
                if file_ext=='.SMILES':
                        name=test
                else:
                        name = file_name
        # edit lines to match limits and specifications of queuing system
        pbs_lines='#PBS -l walltime=120:00:00'
        pbs_lines=pbs_lines + '#PBS -l nodes='+num_procs+ '\n'
        pbs_lines=pbs_lines + '#PBS -N ' + name + '\n'
        pbs_lines=pbs_lines + '#PBS -o ' + name +'.out\n'
        pbs_lines=pbs_lines + '#PBS -j oe\n\n'
        pbs_lines=pbs_lines + 'cd ' + location + '\n\n'
        # edit lines to 
        pbs_lines=pbs_lines + 'module load openbabel\n'
        pbs_lines=pbs_lines + 'module load intel64/14.0.0.080\n'                ##Version dependant on system##
        pbs_lines=pbs_lines + 'module load openmpi_intel64/1.6.5_intel14\n'     ##Version dependant on system##
        pbs_lines=pbs_lines + 'module load amber\n'
        pbs_lines=pbs_lines + 'module load vmd64\n'
        pbs_lines=pbs_lines + 'module load namd/2.10\n\n'
        pbs_lines=pbs_lines + 'python flux.py -n $[PBS_NUM_NODES*24] -c ' + num_clusters + ' -i "' + file_ext[1:] + ' ' + input+' '+name+ '\n'
        return(pbs_lines)

sublines=genSubScript()

submit = subprocess.Popen(["qsub"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
sublines=PBS_directives+attr+joblines
child_input.write(sublines)
child_input.close()
jobid=child_stdout.read()

submitout, submiterr = submit.communicate()
print("Job submitted to queue as: "+ child_stdout+"\n")
