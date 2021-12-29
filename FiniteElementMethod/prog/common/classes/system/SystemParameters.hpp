#pragma once

/**
* Data class for storing parameters of the solving equation
 * */
class SystemParameters
{
public:
    //number of nodes in each direction
    int n;
    
    //length of domain
    double L;

    //speed of the liquid
    double U;

    //incline of the surface
    double k;

    //viscosity of the liquide
    double mu;

    //beginning of the coordinate system
    double xOrigin;
    double yOrigin;

    //border conditions on values of pressure
    double LOW_BORDER;
    double HIGH_BORDER;

    //characteristic size of pressure for dimensionlessness
    double pMin;


    double borderLength;
    double hMin;
    double Hn;

    SystemParameters()
    {
        n = 0;
        mu = 0.0;
        L = 0.0;
        U = 0.0;
        pMin = 0.0;
        k = 0.0;
        hMin = 0.0;
        HIGH_BORDER = 2.0;
        LOW_BORDER = 1.0;
        xOrigin = 0.0;
        yOrigin = 0.0;
        k = 0.0;
        borderLength = 0.0;
    }
};