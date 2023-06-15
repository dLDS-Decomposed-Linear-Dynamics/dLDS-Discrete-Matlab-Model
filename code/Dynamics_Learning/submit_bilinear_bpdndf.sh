#PBS -N bilinear-bpdndf-test
#PBS -q rozell
#PBS -j oe 
#PBS -m abe
#PBS -l walltime=48:00:00 
#PBS -l nodes=1:ppn=20
#PBS -l mem=50gb 

echo Running on host `uname -n`
echo Time is `date`
echo Directory is `pwd`

cd /nv/hp16/acharles6/VersionedCode/2013_Flearning/code/Dynamics_Learning/
echo Directory is `pwd`

matlab -nodisplay -singleCompThread -r "FORCE_learn_all_bilinear_bpdndf" -logfile /nv/hp16/acharles6/VersionedCode/2013_Flearning/results/FORCE_test.txt 



