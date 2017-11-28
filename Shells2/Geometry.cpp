#include "Geometry.h"

using namespace Eigen;

Matrix2m firstFundamentalForm(const Vector3m &p0, const Vector3m &p1, const Vector3m &p2)
{
    Matrix2m result;
    result << (p1 - p0).dot(p1 - p0), (p1 - p0).dot(p2 - p0),
        (p1 - p0).dot(p2 - p0), (p2 - p0).dot(p2 - p0);
    return result;
}

Matrix3m crossMatrix(const Vector3m &v)
{
    Matrix3m result;
    result << 0, -v[2], v[1],
        v[2], 0, -v[0],
        -v[1], v[0], 0;
    return result;
}