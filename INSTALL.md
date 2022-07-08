## To compile and install the file

Step 1: rememeber to change the following path

1.1. PYTHON_DIR, the python that you have installed
1.2. NUMPY_DIR, the corresponding numpy library that you have installed, we need the dir that you can find the numpy header files
1.3. INSTALL_DIR, the location where you wanted that library to be installed, normally to find the anaconda env location

### COMPILING DEPENDENCY
#### PYTHON_DIR=/home/hongfeng/anaconda3/envs/thinfilm/include/python3.8/
#### NUMPY_DIR=/home/hongfeng/anaconda3/envs/thinfilm/lib/python3.8/site-packages/numpy/core/include/

### THE PYTHON LIBRARY YOU WANT TO BE INSTALLED
#### INSTALL_DIR=/home/hongfeng/anaconda3/envs/thinfilm/lib/python3.8/site-packages/

Step 2: compiling and install

it is then very straightforward to do

### make dependency multilater.o
#### 2.1. cd to where you can find multilayer.cpp and multilayer.h, then "make" 
#### 2.2. cd to python where you can find pymultilayer.i, then "make"
