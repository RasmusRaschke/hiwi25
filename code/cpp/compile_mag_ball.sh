# shh to the compile node and clone from github.
ssh compile22-amd-1
module load git
cd /superscratch/rraschke/
git clone https://github.com/RasmusRaschke/hiwi25.git
cd hiwi25/cpp/src

# load petsc for eigen (if you use a intel node you need to load a different module)
module purge
module load petsc/3.21.1-real-4amd
g++ -o sim.out -I/usr/include/eigen3 sim.cpp
./sim.out < input 