ddtFoam, an OpenFOAM Solver to simulate the deflagration-to-detonation transition in gas mixtures

For an introduction to the solver and tutorials, see the file readme.pdf

**MATLAB files readme**

2016-03-28

Added the files lookupY_05.m and lookuptIgn_06.m.
These are (Matlab based) Cantera files which you can use to generate your own mixture composition and ignition delay times.
They require the installation of the Shock and Detonation Toolbox, see
http://shepherd.caltech.edu/EDL/public/cantera/html/SD_Toolbox/

**Build**

Find installed OpenFoam `bashrc` by running `find / -name bashrc` command, Then run: 

```
source OpenFOAM_PATH_TO/etc/bashrc
```

Then run this command in the cloned repository directory:

```
./Allwmake
```

**Build in minimum OpenFoam 2.1.1 docker**

Follow instructions in https://github.com/rezamarzban/openfoam211-docker to install and run its docker.

Then run below commands in the `~/OpenFOAM/root-2.1.x` directory:

```
git clone https://git.code.sf.net/p/ddtfoam/code ddtfoam-code
cd ddtfoam-code
./Allwmake
```

**Simulation example in the docker**

Run:

```
cd /OpenFOAM/root-2.1.x/ddtFoam-code/tutorials/pddtFoam_Tutorial
./Setup
pddtFoam
```
