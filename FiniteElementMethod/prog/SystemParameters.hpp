/**
 * This class stores input data
 * */
class SystemPatemeters
{
public:
    int n;
    double L;
    double U;
    double k;
    double mu;
    double Hn;
    double hMin;
    double pMin;
    double xOrigin;
    double yOrigin;
    double LOW_BORDER;
    double HIGH_BORDER;
    double borderLength;

    SystemPatemeters()
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