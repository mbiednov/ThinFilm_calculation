#ifndef multilayer_h
#define multilayer_h

#include <math.h>
#include <iostream>
#include <fstream>

const double Pi= 3.1415926535897932385;
const double condVac=1/377.0;

class mycomplex {
    
public:
    
    double Re,Im;
    mycomplex();
    mycomplex(double a, double b);
    
    void setvalue(double a, double b);
    void setvalue(mycomplex number);
    
    mycomplex operator + (mycomplex number);
    mycomplex operator + (double number);
    mycomplex operator * (mycomplex number);
    mycomplex operator * (double number);
    mycomplex operator / (mycomplex number);
    mycomplex operator / (double number);
    mycomplex operator - (mycomplex number);
    mycomplex operator - (double number);
    mycomplex operator = (mycomplex number);
};

struct Ellipsommetry
{
    //Structure to hold the ellipsommetry values tg(Psi) and cos(Delta)
    double tgPsi;
    double cosDelta;
};

struct Period
{
    //This structure will contain the information about a Period
    //A period is a sequence of optical layers that enters the multilayer strcture
    //more than one time in a row
    int numberOfRepetitions; //number of times that this sequence is included in the system
    mycomplex M[2][2]; //Matrix of the period

};

double absval(mycomplex& number);
double absvalsq(mycomplex& number);
mycomplex  complexsqrt(mycomplex number);
mycomplex complexcos(mycomplex z);
mycomplex complexsin(mycomplex z);
mycomplex imultiplication(mycomplex number);
mycomplex complexconj(mycomplex number);

void printComplexMatrix(mycomplex M[2][2]);
void printComplexVector(mycomplex V[2]);

//***********
//functions needed to build matrix of the layer, matrix of the entire system and so on
//these functions are called from "basic" or "main" calculation functions(next section of functions)


//multiplication of two complex matrices of the size 2x2. The result is saved to the result matrix
void multiplyComplexMatr(mycomplex matr1[2][2], mycomplex matr2[2][2], mycomplex result[2][2]);

//Nultiplication of the complex matrix 2x2 by the complex vector
void multiplyComplexMatrVector(mycomplex M[2][2], mycomplex V[2], mycomplex result[2]);

//Creation of the matrix of one layer of the system with given parameters, at given wavelength.
//Complex sin and cos of the angle of the beam inside the layer have to be provided to the function
void makeLayersMatrix(mycomplex n, double thickness, double wavelength, mycomplex layersCosine, mycomplex layersSine, bool pPolarized, mycomplex layersMatrix[2][2]);

//Creation of the complex vector that represents the substrate
void makeSubstrateVector(int numberOfLayers , double n[], double k[],double wavelength, double angleOfIncidence, bool pPolarized, mycomplex result[2]);

//Creation of the matrix of the layer system as a product of matriced of all layers
//angleOfIncidence is the AOI at the first interface n[0]/n[1]
//k[0] has to be 0, h of the first and last layers is 0 too.
void makeMatrixOfTheLayerSystem(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized, mycomplex finalResult[2][2]);

//Creation of the matrix of the layer system if it specifyed as a sequence of periods.
//Same as previous function
//The first layer of the total system (medium) can be represented as a period that consists of one layer with one repetition
void makeMatrixOfTheLayerSystem(Period periods[], int numberOfPeriods, mycomplex result[2][2]);

//Creation of the period that consists of given layers, repeated numberOfRepetitions times.
//In total, a period has numberOfLayers*numberOfRepetitions layers
//nmedium is the non-complex refraction index of the medium (first layer of the system in previous functions)
//nmedium is needed in order to know the AOI on the first layer of the period
Period createPeriod(int numberOfLayers, int numberOfRepetitions, double nmedium, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized);


//***********
//Basic calculation functions
//these functions are called from "main" calculation functions from the following sections

//Calculation of the complex refletivity of the multilayer system at given AOI, wavelength and polarization
//The multilayer system is specifyed as a squence of the layers with given n,k and thickness
mycomplex calculateComplexReflectivity(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized);

//Calculation of the complex refletivity of the multilayer system at given AOI, wavelength and polarization
//The multilayer system is specifyed as a array of periods. The total number of elements in the vector periods[] must be numberOfPeriods
mycomplex calculateComplexReflectivity(Period periods[], int numberOfPeriods, double nmedium, double nsubstrate, double ksubstrate, double angleOfIncidence, double wavelength, bool pPolarized);


//***********
//Main reflectivty calculation functions
//The function are supprosed to be directly used by the user


//Calculation of the refletivity of the multilayer system at given AOI, wavelength and polarization
//The multilayer system is specifyed as a squence of the layers with given n,k and thickness
double calculateReflectivity(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized);


//Calculation of the refletivity of the multilayer system at given AOI, wavelength and polarization
//The multilayer system is specifyed as a squence of periods
//The periods[] vector has to be created before
double calculateReflectivity(Period periods[], int numberOfPeriods, double nmedium, double nsubstrate, double ksubstrate, double angleOfIncidence, double wavelength, bool pPolarized);

//Calculation of the transmission
double calculateTransmission(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized);

//Calcualtion of the absorption
double calculateAbsorption(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized);


//Calcaultion of the ellipsometric parameters tan(psi) and cos(delta)

//A multiplayer system is specifyed as a sequence of layers
Ellipsommetry calculateEllispommetricValues(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength);

//A multiplayer system is specifyed as a sequence of periods
Ellipsommetry calculateEllispommetricValues(Period periods[], int numberOfPeriods, double nmedium, double nsubstrate, double ksubstrate, double angleOfIncidence, double wavelength, bool pPolarized);




#endif /* multilayer_h */
