/**
 * This class stores input data
 * */
class SystemPatemeters
{
public:
    /**
    * Number of finite elements
    * */
    int n;

    /**
     * Fluid speed
     * */
    double U;

    /**
     * Fluid viscosity
     * */
    double mu;

    /**
     * Upper surface slope parameter
     * */
    double k;

    /**
     * Real coefficient k
     * */
    double realK;

    /**
    * Height of input gap
    * */
    double hMin;

    /**
     * Real value of input gap
     * */
    double realInputGap;

    /**
    * Real pressure on borders
    * */
    double realPressure;

    /**
     * Length of the region
     * */
    double L;

    /**
     * X and Y coordinates of top-left point of the region
     * */
    double xOrigin;
    double yOrigin;

    /**
     * Value of pressure on three sides of the region
     * */
    double LOW_BORDER;

    /**
     * Value of pressure on one side of the region
     * */
    double HIGH_BORDER;

    SystemPatemeters()
    {
        n = 0;
        U = 0.0;
        L = 0.0;
        mu = 0.0;
        xOrigin = 0.0;
        yOrigin = 0.0;

        k = 0.0;
        hMin = 0.0;
        HIGH_BORDER = 2.0;
        LOW_BORDER = 1.0;
        
        realK = 0.0;
        realInputGap = 0.0;
        realPressure = 0.0;
    }
};