#ifndef ARBITRARYPRECISION_H
#define ARBITRARYPRECISION_H

#include <Eigen/Core>
#include "unsupported/Eigen/MPRealSupport"


typedef mpfr::mpreal scalar;

//typedef double scalar;

typedef Eigen::Matrix<scalar, 3, 1> Vector3m;
typedef Eigen::Matrix<scalar, 2, 1> Vector2m;
typedef Eigen::Matrix<scalar, 2, 2> Matrix2m;
typedef Eigen::Matrix<scalar, 3, 3> Matrix3m;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> VectorXm;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 3> MatrixX3m;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXm;

#endif
