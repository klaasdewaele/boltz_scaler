!/bin/bash  
#PBS -N boltz_prediction_{complex_num}  
#PBS -l nodes=1:ppn=8,gpus=1  
#PBS -l walltime=0:30:0  
#PBS -l mem=64gb  

module load Boltz-1/0.4.1-foss-2023a-CUDA-12.1.1  
