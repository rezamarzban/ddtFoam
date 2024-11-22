ddtFoam, an OpenFOAM Solver to simulate the deflagration-to-detonation transition in gas mixtures

For an introduction to the solver and tutorials, see the file readme.pdf

**MATLAB files readme**

2016-03-28

Added the files lookupY_05.m and lookuptIgn_06.m.
These are (Matlab based) Cantera files which you can use to generate your own mixture composition and ignition delay times.
They require the installation of the Shock and Detonation Toolbox, see
http://shepherd.caltech.edu/EDL/public/cantera/html/SD_Toolbox/

It needs `OpenFoam 2.1.1` with an old OS version to building.

**Run with Docker**

https://hub.docker.com/r/rm314159/ddtfoam

Pull: 

```
docker pull rm314159/ddtfoam:latest
```

Install: 

```
docker run  -it -d --name container_ddt --entrypoint /bin/bash rm314159/ddtfoam
```

Run: 

```
docker start -i container_ddt
```

Everytime in the running container, Run:

```
source /opt/openfoam211/etc/bashrc
```

Simulation example, Run:

```
cd ~/OpenFOAM/root-2.1.x/ddtFoam-code/tutorials/pddtFoam_Tutorial
./Setup
pddtFoam
```

**Build with Docker**

To install and run its docker image and container, Follow instructions in:

https://github.com/rezamarzban/openfoam211-docker

or:

https://github.com/rezamarzban/openfoam211

Then run below commands in the `~/OpenFOAM/root-2.1.x` directory:

```
git clone https://git.code.sf.net/p/ddtfoam/code ddtfoam-code
cd ddtfoam-code
./Allwmake
```

**Recommendation**

While `OpenFoam` and its solvers is not user-friendly, Use `CONVERGE Studio CFD` and its solvers which is user-friendly and available for different `Microsoft Windows` OS to simulating DDT. Also there is `ANSYS FLUENT`, `COMSOL Multiphysics` and ` Autodesk CFD` to simulating such this problem.
