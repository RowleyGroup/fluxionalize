import sys, getopt, os, subprocess

#write test.sh and send it to queue
def write_sh():
        location = os.getcwd()
        name_spec=False
        num_procs='2'
        num_clusters='1'
        opts, args = getopt.getopt(sys.argv[1:], "p:n:l:i:c:")
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

                        file_name, file_ext=os.path.splitext(input)

                        if not any (file_ext in type for type in file_formats):
                                print("File type %s may not be supported.") %file_ext
                        if file_ext=="":
                                file_ext='.SMILES'

        if name_spec==False:
                if file_ext=='.SMILES':
                        name=test
                else:
                        name = file_name
        file_out=open(name+'.sh', 'w')
        file_out.write('#PBS -l walltime=120:00:00\n#PBS -l nodes='+num_procs+'\n#PBS -N '+name+'\n#PBS -o '+name+'.out\n#PBS -j oe\n\ncd '+location+'\n\n'
                        +'module load openbabel\nmodule load intel64/14.0.0.080\nmodule load openmpi_intel64/1.6.5_intel14\nmodule load amber\nmodule load vmd64\nmodule load namd/2.10\n\n'
                        +'python tes.py -n $[PBS_NUM_NODES*24] -c '+num_clusters+' -i "'+file_ext[1:]+' '+input+' '+name+'"\n')
        file_out.close()


write_sh()
submit = subprocess.Popen("qsub "+name+".sh", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
submitout, submiterr = submit.communicate()
print("Job submitted to queue as: "+submitout+"\n")
