#include <iostream>

class Pymultilayer
{
public:

    Pymultilayer();
    ~Pymultilayer();
    
    double calReflectivity(int numberOfLayers, double *n, int Nx1, double *k, int Nx2, double *thickness, int Nx3, double angleOfIncidence, double wavelength, int pPolarized);

    double calTransmission(int numberOfLayers, double *n, int Nx1, double *k, int Nx2, double *thickness, int Nx3, double angleOfIncidence, double wavelength, bool pPolarized);

private:

    
};
