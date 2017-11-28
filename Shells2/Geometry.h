#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "ArbitraryPrecision.h"

Matrix2m firstFundamentalForm(const Vector3m &p0, const Vector3m &p1, const Vector3m &p2);
Matrix3m crossMatrix(const Vector3m &v);

#endif
