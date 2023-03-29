# ThinFilm_calculation
C++ code to perform calculation of the optical properties of thin films.

## Motivation
The goal is to create a C++ library that is capable of calculating reflection, transmission, absorption and ellipsometric 
parameters (tg(Psi) and cos(Delta)) of any given multilayer thin film stack, including stacks that consist of periodically repeating layers.

## Applicability
The code can be used by people who design various types of optical coatings, like anti-reflection coatings, mirrors, beamsplitters, etc. 
In addition it may be usefull to
those, who analyze/simulate ellipsometric data.

## Current structure
The projects contains two files:
* multilayer.h
* multilayer.cpp

## Building the project
* clone the repository to yor own project
* include **multilayer.h** in you project

## Examples of usage
More information on the code and examples of its usage can be found in the **Manual.pdf** in the **docs** folder.  
In addition, the code (without the support of periodically repeating layers) is implemented in the Opal software, that can be found under https://github.com/mbiednov/opal/releases

### Python example: a 100nm thin layer (index 2.22) sandwitched with air (index 1) and glass (index 1.52)
#### Compile:
1. install "swig" for python.
2. open windows cmd or linux terminal and enter into the "python" folder, "swig -c++ -python pymultilayer.i". it will generate "Pymultilayer_wrap.cxx" and "pymultilayer.py".
3. for windows users: open visual studio and create a new project, copy "Pymultilayer_wrap.cxx", "pymultilayer.cpp", and "pymultilayer.h" into the project. compile them into a dll file. renamed the generated dll file as "\_pymultilayer.pyd". 
4. for windows users: copy files of "pymultilayer.py" and "\_pymultilayer.pyd" into a new folder, these two files is a complete library.

#### compile for Linux
1. type "make" at the folder containing "multilayer.cpp", it would generate dependent lib into python folder
2. enter python folder, change 
```
# COMPILING DEPENDENCY
PYTHON_DIR=/home/username/anaconda3/envs/thinfilm/include/python3.8/ # CHANGE to your python path, I use conda
NUMPY_DIR=/home/username/anaconda3/envs/thinfilm/lib/python3.8/site-packages/numpy/core/include/ # CHANGE to your numpy path

# THE PYTHON LIBRARY YOU WANT TO BE INSTALLED
INSTALL_DIR=/home/username/anaconda3/envs/thinfilm/lib/python3.8/site-packages/ # CHANGE to your python site-packages path
```

The results are compared with FDTD. Due to the increased time in simulating reflectity by FDTD, the example only compares the range of incidence angle from 0 to 40 deg. The wavelength is 632nm.
<p align="center">
<img src="https://github.com/MarkMa1990/ThinFilm_calculation/blob/working/docs/comparison.png" width="70%" height="70%">
</p>

### Matlab example:
#### compile: 
1. setup mex compile environment "mex -setup". if you don't know how to setup, just goole it.
2. enter "matlab" folder containing the "mexFilm.cpp" file, type "mex mexFilm.cpp" in matlab cmd. 
3. in matlab cmd, "test" and press "Enter" key to run.

The results are compared with Essential Macleod. 

<p align="center">
<img src="https://github.com/MarkMa1990/ThinFilm_calculation/blob/working/matlab/matlabFigure.jpg" width="70%" height="70%">
</p>
