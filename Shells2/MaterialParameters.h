#ifndef MATERIALPARAMETERS_H
#define MATERIALPARAMETERS_H

struct MaterialParameters
{
    MaterialParameters() : YoungsModulus(1.0e8), PoissonsRatio(0.3) {}

    double YoungsModulus;
    double PoissonsRatio;

    double LameAlpha() const
    {
        return YoungsModulus * PoissonsRatio / (1.0 + PoissonsRatio) / (1.0 - 2.0*PoissonsRatio);
    }

    double LameBeta() const
    {
        return YoungsModulus / (2.0 * (1.0 + PoissonsRatio));
    }
};

#endif