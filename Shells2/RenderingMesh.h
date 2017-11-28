#ifndef RENDERINGMESH_H
#define RENDERINGMHES_H

#include "Eigen/Core"
#ifdef GUI_VERSION
#include <mutex>
#endif
class RenderingMesh
{
public:
    RenderingMesh();
    void recomputeColors();
    void writeMesh(int number);

    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
    Eigen::MatrixX3d C;

    bool showGradientDir;
    Eigen::VectorXd gradientDir;
    bool showDescentDir;
    Eigen::VectorXd descentDir;
    double gradScale;
    Eigen::VectorXd faceVals;
    double colorMax;

    std::string outputPath;
#ifdef GUI_VERSION
    std::mutex visMutex; 
#endif   
};

#endif
