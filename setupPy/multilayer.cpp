#include <stdio.h>
#include "multilayer.h"

mycomplex::mycomplex(){Re=0; Im=0;}
mycomplex::mycomplex(double a, double b){Re=a; Im=b;}
mycomplex mycomplex::operator + (mycomplex number)
{
    mycomplex result;
    result.Re = Re + number.Re;
    result.Im = Im + number.Im;
    return result;
}
mycomplex mycomplex::operator + (double number)
{
    mycomplex result;
    result.Re = Re + number;
    result.Im = Im;
    return result;
}
mycomplex mycomplex::operator * (mycomplex number)
{
    mycomplex result;
    result.Re = Re*number.Re - Im*number.Im;
    result.Im = Re*number.Im + Im*number.Re;
    return result;
}
mycomplex mycomplex::operator * (double number)
{
    mycomplex result;
    result.Re = Re*number;
    result.Im = Im*number;
    return result;
}
mycomplex mycomplex::operator / (mycomplex number)
{
    mycomplex result;
    double r = number.Re*number.Re + number.Im*number.Im;
    result.Re = (Re*number.Re + Im*number.Im) / r;
    result.Im = (Im*number.Re - Re*number.Im) / r;
    return result;
}

mycomplex mycomplex::operator / (double number)
{
    mycomplex result;
    result.Re = Re / number;
    result.Im = Im / number;
    return result;
}
mycomplex mycomplex::operator - (mycomplex number)
{
    mycomplex result;
    result.Re = Re - number.Re;
    result.Im = Im - number.Im;
    return result;
}
mycomplex mycomplex::operator - (double number)
{
    mycomplex result;
    result.Re = Re - number;
    result.Im = Im;
    return result;
}
mycomplex mycomplex::operator = (mycomplex number)
{
    Re = number.Re;
    Im = number.Im;
    return *this;
}
void mycomplex::setvalue (double a, double b)
{
    Re = a;
    Im = b;
}
void mycomplex::setvalue(mycomplex number)
{
    Re = number.Re;
    Im = number.Im;
}



double absval(mycomplex& number)
{
    return sqrt(number.Re*number.Re + number.Im*number.Im);
}

double absvalsq(mycomplex& number)
{
    return (number.Re*number.Re + number.Im*number.Im);
}

mycomplex  complexsqrt(mycomplex number)
{
    double m, f;
    mycomplex result;
    m = sqrt(number.Re*number.Re + number.Im*number.Im);
    f = atan2(number.Im,number.Re);
    result.Re = sqrt(m)*cos(f / 2);
    result.Im = sqrt(m)*sin(f / 2);
    return result;
}

mycomplex complexcos(mycomplex z)
{
    mycomplex result;
    result.Re = cos(z.Re)*(exp(-z.Im) + exp(z.Im)) / 2;
    result.Im = sin(z.Re)*(exp(-z.Im) - exp(z.Im)) / 2;
    return result;
}

mycomplex complexsin(mycomplex z)
{
    mycomplex result;
    result.Re = sin(z.Re)*(exp(-z.Im) + exp(z.Im)) / 2;
    result.Im = -cos(z.Re)*(exp(-z.Im) - exp(z.Im)) / 2;
    return result;
}

mycomplex imultiplication(mycomplex number)
{
    mycomplex result;
    result.Re = -number.Im;
    result.Im = number.Re;
    return result;
}

mycomplex complexconj(mycomplex number)
{
    mycomplex result;
    result.Re=number.Re;
    result.Im=-number.Im;
    return result;
}

void printComplexMatrix(mycomplex M[2][2])
{
    for(int row=0; row<2; row++)
    {
        for (int column=0; column<2; column++)
        {
            std::cout<<M[row][column].Re<<"+i*"<<M[row][column].Im<<"       ";
        }
        std::cout<<std::endl;
    }
}

void printComplexVector(mycomplex V[2])
{
   std::cout<<V[0].Re<<"+i*"<<V[0].Im<<"   ";
    std::cout<<V[1].Re<<"+i*"<<V[1].Im<<std::endl;
}



void multiplyComplexMatr(mycomplex a[2][2], mycomplex b[2][2], mycomplex result[2][2])
{
    //this function multiplies two matrices a*b and saves the result in result
    result[0][0].setvalue(a[0][0]*b[0][0]+a[0][1]*b[1][0]);
    result[0][1].setvalue(a[0][0]*b[0][1]+a[0][1]*b[1][1]);
    result[1][0].setvalue(a[1][0]*b[0][0]+a[1][1]*b[1][0]);
    result[1][1].setvalue(a[1][0]*b[0][1]+a[1][1]*b[1][1]);
}

void multiplyComplexMatrVector(mycomplex M[2][2], mycomplex V[2], mycomplex result[2])
{
    //this function multiplyes a compelx matrix[2][2] by a complex vector and saves the result in the
    //result[2]
    result[0].setvalue(M[0][0]*V[0]+M[0][1]*V[1]);
    result[1].setvalue(M[1][0]*V[0]+M[1][1]*V[1]);
}


void makeLayersMatrix(mycomplex n, double thickness, double wavelength, mycomplex layersCosine, mycomplex layersSine, bool pPolarized, mycomplex layersMatrix[2][2])
{
    //This function fills layersMatrix wih the values according to the optical properties of the layer
    //and its thickness
    //pPolarized shows whether the beam is p-polarized
    //if not, then s-polarization is assumed
    
    mycomplex cond(0,0);
    mycomplex delta(0,0);
    
    if (pPolarized)
    {
        //p-polarization
        cond=n*condVac/layersCosine;
    }
    else
    {
        //s-polarization
        cond=n*condVac*layersCosine;
    }
    
    delta=n*layersCosine*2*Pi*thickness/wavelength;
    
    //setting the values of the matrix elements
    layersMatrix[0][0].setvalue(complexcos(delta));
    layersMatrix[0][1].setvalue(imultiplication(complexsin(delta))/cond);
    layersMatrix[1][0].setvalue(imultiplication(cond*complexsin(delta)));
    layersMatrix[1][1].setvalue(complexcos(delta));
    
}

void makeSubstrateVector(int numberOfLayers , double n[], double k[],double wavelength, double angleOfIncidence, bool pPolarized, mycomplex result[2])
{
    //This function makes the vector of the substrate and saves it to the result[2]
    //the result has a form [1 cond_substrate]
    mycomplex sinAngle(0,0);
    mycomplex cosAngle(0,0);
    mycomplex cond_substrate(0,0);
    
    mycomplex N0(n[0],k[0]);
    mycomplex Ns(n[numberOfLayers-1],k[numberOfLayers-1]);
    
    //calculating the angle of the beam inside the substrate
    sinAngle=N0*sin(angleOfIncidence)/Ns;
    cosAngle=complexsqrt(sinAngle*sinAngle*(-1) + 1 );
    
    if (pPolarized)
    {
        //p-polarization
        cond_substrate=Ns*condVac/cosAngle;
    }
    else
    {
        //s-polarization
        cond_substrate=Ns*condVac*cosAngle;
    }
    result[0].setvalue(1, 0);
    result[1].setvalue(cond_substrate);
}

void makeMatrixOfTheLayerSystem(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized, mycomplex finalResult[2][2])
{
    //This function calculates the rmatrix of a multilayer system
    //angle of incidence and  at given wavelength.
    //n[] contains values for the real part of refraction index of each layer
    //k[] contains the values for the imaginary part of the refraction index correspondingly to n[]
    //pPolarized defines whether the light is p-polarized, else s-polarization is considered
    
    mycomplex N(0,0);
    mycomplex N0(n[0],k[0]);
    mycomplex sinAngle(0,0); //values for the complex sine and cosine of the angle of incidence in the
    mycomplex cosAngle(0,0); //current layer
    mycomplex intermediateResult[2][2]; //matrix to hold the intermediate calculation values
    mycomplex currentLayerMatr[2][2]; //matrix for the current layer of the calculation
    
    //initailizing matrix elements
    for(int row=0; row<2; row++)
    {
        for (int column=0; column<2; column++)
        {
            intermediateResult[row][column].setvalue(0, 0);
            if (row==column)
            {
                finalResult[row][column].setvalue(1, 0);
            }
            else
            {
                finalResult[row][column].setvalue(0, 0);
            }
        }
    }
    
    //now we start the calculation loop
    //within the loop matrix of the system is calculated step by step by multiplying finalResult by the
    //currentLayer from the right side
    for (int currentLayer=0; currentLayer<numberOfLayers; currentLayer++)
    {
        //We iterate through all matrices though in the calculation only matrices of layers except for the
        //first one (medium) and the last one (substrate) are needed. Due to the fact that the thicknes
        //of these two layers is 0, they are the unity matrices and therefore don't affect the result
        //first we create matrix of the current layer
        N.setvalue(n[currentLayer], k[currentLayer]); //complex refraction index of the current layer
        sinAngle=N0*sin(angleOfIncidence)/N;
        cosAngle=complexsqrt(sinAngle*sinAngle*(-1) + 1 );
        
        makeLayersMatrix(N, thickness[currentLayer], wavelength, cosAngle, sinAngle, pPolarized, currentLayerMatr);
        
        multiplyComplexMatr(finalResult, currentLayerMatr, intermediateResult);
        //now matrix of the processed layers is intermediateResult
        //we save intermediate result to the final result
        for (int row=0; row<2; row++)
        {
            for (int column=0; column<2; column++)
            {
                finalResult[row][column].setvalue(intermediateResult[row][column]);
            }
        }
        
    }

    
    
}


//Calculation functions

mycomplex calculateComplexReflectivity(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized)
{
    //This functions calculates complex reflectivity of a given multilayer stack at a given angle of incidence
    //and at a given wavelength for a given polarization
    
    //Converting the AOI to radians;
    double angleOfIncidenceRad=angleOfIncidence*Pi/180;
    mycomplex result(0,0);
    mycomplex Y(0,0);
    mycomplex r(0,0); //compelx reflectivity
    mycomplex matrixOfTheSystem[2][2];
    mycomplex vectorOfTheSubstrate[2];
    mycomplex vectorOfTheSystem[2]; //Vector the holds the result of multplication of M of the System by the
    //Vector of the substrate
    double conductivityMedium; //conductivity of the medium (first layer)
    
    //initializing the matrix
    for (int row=0; row<2; row++)
    {
        for (int column=0; column<2; column++)
        {
            matrixOfTheSystem[row][column].setvalue(0, 0);
        }
    }
    //initializing the vector
    vectorOfTheSubstrate[0].setvalue(0, 0);
    vectorOfTheSubstrate[1].setvalue(0, 0);
    vectorOfTheSystem[0].setvalue(0, 0);
    vectorOfTheSystem[1].setvalue(0, 0);
    
    //Starting the calculation
    makeMatrixOfTheLayerSystem(numberOfLayers, n, k, thickness, angleOfIncidenceRad, wavelength, pPolarized, matrixOfTheSystem);
    makeSubstrateVector(numberOfLayers, n, k, wavelength, angleOfIncidenceRad, pPolarized, vectorOfTheSubstrate);
    
    multiplyComplexMatrVector(matrixOfTheSystem, vectorOfTheSubstrate, vectorOfTheSystem);
    
    //conductivity of the medium
    if (pPolarized)
    {
        //p-polarization
        conductivityMedium=n[0]*condVac/cos(angleOfIncidenceRad);
    }
    else
    {
        //s-polarization
        conductivityMedium=n[0]*condVac*cos(angleOfIncidenceRad);
    }
    
    Y.setvalue(vectorOfTheSystem[1]/vectorOfTheSystem[0]);
    result=(Y*(-1)+conductivityMedium)/(Y+conductivityMedium);
    
    
    return result;
}

mycomplex calculateComplexReflectivity(Period periods[], int numberOfPeriods, double nmedium, double nsubstrate, double ksubstrate, double angleOfIncidence, double wavelength, bool pPolarized)
{
    //This function calculates complex reflectivity of the multilayer system, stored in periods[]
    //In periods, each element is a period structure with two fields: number of repetitions and the matrix of the
    //period
    mycomplex result(0,0);
    mycomplex Y(0,0);
    mycomplex sinAngle(0,0);
    mycomplex cosAngle(0,0);
    mycomplex N0(nmedium,0);
    mycomplex Ns(nsubstrate,ksubstrate);
    mycomplex cond_substrate(0,0);
    mycomplex cond_medium(0,0);
    mycomplex matrixOfTheSystem[2][2];
    mycomplex vectorOfTheSubstrate[2];
    mycomplex vectorOfTheSystem[2];
    //Initializing the vectors
    vectorOfTheSubstrate[0].setvalue(0, 0);
    vectorOfTheSubstrate[1].setvalue(0, 0);
    vectorOfTheSystem[0].setvalue(0, 0);
    vectorOfTheSystem[1].setvalue(0, 0);
    
    
    sinAngle=N0*sin(angleOfIncidence*Pi/180)/Ns;
    cosAngle=complexsqrt(sinAngle*sinAngle*(-1) + 1 );
    
    if (pPolarized)
    {
        //p-polarization
        cond_substrate=Ns*condVac/cosAngle;
    }
    else
    {
        //s-polarization
        cond_substrate=Ns*condVac*cosAngle;
    }
    vectorOfTheSubstrate[0].setvalue(1, 0);
    vectorOfTheSubstrate[1].setvalue(cond_substrate);
    
    makeMatrixOfTheLayerSystem(periods, numberOfPeriods, matrixOfTheSystem);
    
    multiplyComplexMatrVector(matrixOfTheSystem, vectorOfTheSubstrate, vectorOfTheSystem);
    
    //conductivity of the medium
    if (pPolarized)
    {
        //p-polarization
        cond_medium=N0*condVac/cos(angleOfIncidence*Pi/180);
    }
    else
    {
        //s-polarization
        cond_medium=N0*condVac*cos(angleOfIncidence*Pi/180);
    }
    
    Y.setvalue(vectorOfTheSystem[1]/vectorOfTheSystem[0]);
    result=(Y*(-1)+cond_medium)/(Y+cond_medium);
    
    
    return result;
}

void makeMatrixOfTheLayerSystem(Period periods[], int numberOfPeriods, mycomplex result[2][2])
{
    //This function makes a total matrix of the system that consists of periods[]
    
    //Intializing the result and intermResult
    result[0][0].setvalue(1, 0);
    result[0][1].setvalue(0, 0);
    result[1][0].setvalue(0, 0);
    result[1][1].setvalue(1, 0);
    
    mycomplex intermResult[2][2];
    for (int row=0; row<2; row++)
    {
        for (int column=0; column<2; column++)
        {
            intermResult[row][column].setvalue(0, 0);
        }
    }
    
    
    for (int period=0; period<numberOfPeriods; period++)
    {
        for (int repetition=0; repetition<periods[period].numberOfRepetitions; repetition++)
        {
            multiplyComplexMatr(result, periods[period].M, intermResult);
            //Now we save the intermResult to the result
            for (int row=0; row<2; row++)
            {
                for (int column=0; column<2; column++)
                {
                    result[row][column].setvalue(intermResult[row][column]);
                }
            }
        }
    }
}

Period createPeriod(int numberOfLayers, int numberOfRepetitions, double nmedium, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized)
{
    //This function creates a matrix of the period which is a product of matrices of the layers of this period
    //here angle of incidence is in degrees
    //number of layers is the number of layers, that are being repeatedly added to the system
    Period result;
    mycomplex matrixOfThePeriod[2][2];
    double angleOfIncidenceRad=angleOfIncidence*Pi/180;
    //now we create a new vector for n[],k[] and h[] to modify them, to able to call makeMatrixOfTheLayerSystem
    double *localN, *localK, *localH;
    localN = new double [numberOfLayers+2];
    localK = new double [numberOfLayers+2];
    localH = new double [numberOfLayers+2];
    //first element is the medium, last one is the substrate
    localN[0]=nmedium;
    localN[numberOfLayers+1]=1;
    localK[0]=0;
    localK[numberOfLayers+1]=0;
    localH[0]=0;
    localH[numberOfLayers+1]=0;
    for (int layer=1; layer<=numberOfLayers; layer++)
    {
        localN[layer]=n[layer-1];
        localK[layer]=k[layer-1];
        localH[layer]=thickness[layer-1];
    }
    
    //initializing the matrix
    for (int row=0; row<2; row++)
    {
        for (int column=0; column<2; column++)
        {
            matrixOfThePeriod[row][column].setvalue(0, 0);
        }
    }
    
    result.numberOfRepetitions=numberOfRepetitions;
    makeMatrixOfTheLayerSystem(numberOfLayers+2, localN, localK, localH, angleOfIncidenceRad, wavelength, pPolarized, matrixOfThePeriod);
    //Assigning the values to the result.M
    for(int row=0;row<2; row++)
    {
        for(int column=0; column<2; column++)
        {
            result.M[row][column]=matrixOfThePeriod[row][column];
        }
    }
    
    return result;
}


double calculateReflectivity(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized)
{
    //This function calcualtes the absolute value of a reflectivity
    mycomplex r(0,0);
    double reflectivity{0};
    r=calculateComplexReflectivity(numberOfLayers, n, k, thickness, angleOfIncidence, wavelength, pPolarized);
    reflectivity=absvalsq(r);
    
    return reflectivity;
}




double calculateReflectivity(Period periods[], int numberOfPeriods, double nmedium, double nsubstrate, double ksubstrate, double angleOfIncidence, double wavelength, bool pPolarized)
{
    //This function calculates absolute value of the reflectivity of a periodic structure
    mycomplex r(0,0);
    double reflectivity{0};
    r=calculateComplexReflectivity(periods, numberOfPeriods, nmedium, nsubstrate, ksubstrate, angleOfIncidence, wavelength, pPolarized);
    reflectivity=absvalsq(r);
    
    return reflectivity;
}

//Transmission calculation

double calculateTransmission(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized)
{
    //This functions calculates transmission of a given multilayer stack at a given angle of incidence
    //and at a given wavelength for a given polarization
    
    //Converting the AOI to radians;
    double angleOfIncidenceRad=angleOfIncidence*Pi/180;
    double result{0};
    mycomplex B(0,0);
    mycomplex C(0,0);
    mycomplex Y(0,0);
    mycomplex matrixOfTheSystem[2][2];
    mycomplex vectorOfTheSubstrate[2];
    mycomplex vectorOfTheSystem[2]; //Vector the holds the result of multplication of M of the System by the
    //Vector of the substrate
    double conductivityMedium; //conductivity of the medium (first layer)
    
    //initializing the matrix
    for (int row=0; row<2; row++)
    {
        for (int column=0; column<2; column++)
        {
            matrixOfTheSystem[row][column].setvalue(0, 0);
        }
    }
    //initializing the vector
    vectorOfTheSubstrate[0].setvalue(0, 0);
    vectorOfTheSubstrate[1].setvalue(0, 0);
    vectorOfTheSystem[0].setvalue(0, 0);
    vectorOfTheSystem[1].setvalue(0, 0);
    
    //Starting the calculation
    makeMatrixOfTheLayerSystem(numberOfLayers, n, k, thickness, angleOfIncidenceRad, wavelength, pPolarized, matrixOfTheSystem);
    makeSubstrateVector(numberOfLayers, n, k, wavelength, angleOfIncidenceRad, pPolarized, vectorOfTheSubstrate);
    
    multiplyComplexMatrVector(matrixOfTheSystem, vectorOfTheSubstrate, vectorOfTheSystem);
    
    //conductivity of the medium
    if (pPolarized)
    {
        //p-polarization
        conductivityMedium=n[0]*condVac/cos(angleOfIncidenceRad);
    }
    else
    {
        //s-polarization
        conductivityMedium=n[0]*condVac*cos(angleOfIncidenceRad);
    }
    
    B.setvalue(vectorOfTheSystem[0]);
    C.setvalue(vectorOfTheSystem[1]);
    Y=B*conductivityMedium+C;
    
    result=4*conductivityMedium*vectorOfTheSubstrate[1].Re/absvalsq(Y);
    
    return result;
}

//Absorption calculation

double calculateAbsorption(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength, bool pPolarized)
{
    //This functions calculates Absorption of a given multilayer stack at a given angle of incidence
    //and at a given wavelength for a given polarization
    
    //Converting the AOI to radians;
    double angleOfIncidenceRad=angleOfIncidence*Pi/180;
    double result{0};
    mycomplex B(0,0);
    mycomplex C(0,0);
    mycomplex Y(0,0);
    mycomplex Y2(0,0);
    mycomplex matrixOfTheSystem[2][2];
    mycomplex vectorOfTheSubstrate[2];
    mycomplex vectorOfTheSystem[2]; //Vector the holds the result of multplication of M of the System by the
    //Vector of the substrate
    double conductivityMedium; //conductivity of the medium (first layer)
    
    //initializing the matrix
    for (int row=0; row<2; row++)
    {
        for (int column=0; column<2; column++)
        {
            matrixOfTheSystem[row][column].setvalue(0, 0);
        }
    }
    //initializing the vector
    vectorOfTheSubstrate[0].setvalue(0, 0);
    vectorOfTheSubstrate[1].setvalue(0, 0);
    vectorOfTheSystem[0].setvalue(0, 0);
    vectorOfTheSystem[1].setvalue(0, 0);
    
    //Starting the calculation
    makeMatrixOfTheLayerSystem(numberOfLayers, n, k, thickness, angleOfIncidenceRad, wavelength, pPolarized, matrixOfTheSystem);
    makeSubstrateVector(numberOfLayers, n, k, wavelength, angleOfIncidenceRad, pPolarized, vectorOfTheSubstrate);
    
    multiplyComplexMatrVector(matrixOfTheSystem, vectorOfTheSubstrate, vectorOfTheSystem);
    
    //conductivity of the medium
    if (pPolarized)
    {
        //p-polarization
        conductivityMedium=n[0]*condVac/cos(angleOfIncidenceRad);
    }
    else
    {
        //s-polarization
        conductivityMedium=n[0]*condVac*cos(angleOfIncidenceRad);
    }
    
    B.setvalue(vectorOfTheSystem[0]);
    C.setvalue(vectorOfTheSystem[1]);
    Y=B*conductivityMedium+C;
    Y2=B*complexconj(C)-vectorOfTheSubstrate[1];
    
    result=4*conductivityMedium*Y2.Re/absvalsq(Y);
    
    return result;
}

//Main ellipsommetry calculation functions

Ellipsommetry calculateEllispommetricValues(int numberOfLayers, double n[], double k[], double thickness[], double angleOfIncidence, double wavelength)
{
    //This function calculates ellipsommetric values for a given multilayer stack
    Ellipsommetry result;
    mycomplex Rp(0,0);
    mycomplex Rs(0,0);
    mycomplex z(0,0);
    
    Rp=calculateComplexReflectivity(numberOfLayers, n, k, thickness, angleOfIncidence, wavelength, true);
    Rs=calculateComplexReflectivity(numberOfLayers, n, k, thickness, angleOfIncidence, wavelength, false);
    z=Rp/Rs;
    
    result.tgPsi=absval(z);
    result.cosDelta=z.Re/absval(z);
    
    return result;
}

Ellipsommetry calculateEllispommetricValues(Period periods[], int numberOfPeriods, double nmedium, double nsubstrate, double ksubstrate, double angleOfIncidence, double wavelength, bool pPolarized)
{
    //This function calculates ellipsommetric values for a given multilayer stack
    Ellipsommetry result;
    mycomplex Rp(0,0);
    mycomplex Rs(0,0);
    mycomplex z(0,0);
    
    Rp=calculateComplexReflectivity(periods, numberOfPeriods, nmedium, nsubstrate, ksubstrate, angleOfIncidence, wavelength, true);
    Rs=calculateComplexReflectivity(periods, numberOfPeriods, nmedium, nsubstrate, ksubstrate, angleOfIncidence, wavelength, pPolarized);
    z=Rp/Rs;
    
    result.tgPsi=absval(z);
    result.cosDelta=z.Re/absval(z);
    
    return result;
}


