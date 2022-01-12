#python script to run jobs in IMS-server
import os
import glob
import shutil

njobs=5

listdirs=['Modules','Test_modules']
inputdir='Initial_inputs'
execfile=glob.glob('*.py')

root=os.getcwd()
for i in range(1,njobs+1):
    dirname='job-'+str(i)
    os.makedirs(dirname,exist_ok=True)
    os.chdir(dirname)

    with open('run.sh','w') as js:
        js.write('#!/bin/csh \n')
        js.write('#$ -cwd \n')
        js.write('#$ -V -S /bin/bash \n')
        js.write('#$ -N KMC-job-'+str(i)+'\n')
        js.write('#$ -o stdout \n')
        js.write('#$ -e stderr \n')
        js.write('#$ -q all.q \n')
        js.write('#$ -pe x24 24 \n')
        js.write('export OMP_NUM_THREADS=1 \n')
        js.write('python RFKMC_multi.py \n')

    #copy required

    for dir in listdirs:
        origin=os.path.join(root,dir)
        shutil.copytree(origin, './'+dir)

    #input
    inputname='kmc_input_{0}.yml'.format(i)
    origin=os.path.join(root,inputdir,inputname)
    shutil.copy(origin,'./kmc_input.yml')

    #exec
    for file in execfile:
        origin=os.path.join(root,file)
        shutil.copy(origin,file)

    os.system('qsub run.sh')

    os.chdir(root)
