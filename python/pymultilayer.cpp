#include "pymultilayer.h"

#include "multilayer.h"

Pymultilayer::Pymultilayer()
{

}

Pymultilayer::~Pymultilayer()
{

}

double Pymultilayer::calReflectivity(int numberOfLayers, double *n, int Nx1, double *k, int Nx2, double *thickness, int Nx3, double angleOfIncidence, double wavelength, int pPolarized)
{
    double result=0;
    if (pPolarized==0)
        result = calculateReflectivity(numberOfLayers, n, k, thickness, angleOfIncidence, wavelength, false);
    else
        result = calculateReflectivity(numberOfLayers, n, k, thickness, angleOfIncidence, wavelength, true);
        
    return result;
}

double Pymultilayer::calTransmission(int numberOfLayers, double *n, int Nx1, double *k, int Nx2, double *thickness, int Nx3, double angleOfIncidence, double wavelength, int pPolarized)
{
    double result = 0;
    if(pPolarized==0)
        result = calculateTransmission(numberOfLayers, n, k, thickness, angleOfIncidence, wavelength, false);
    else
        result = calculateTransmission(numberOfLayers, n, k, thickness, angleOfIncidence, wavelength, true);

    return result;

}
