#ifndef SIMULATIONMESH_H
#define SIMULATIONMESH_H

#include "MaterialParameters.h"
#include <Eigen/Core>
#include <vector>
#include "optimization.h"
#include "ArbitraryPrecision.h"

class RenderingMesh;

class SimulationMesh
{
public:
    SimulationMesh(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, double thickness);
    friend SimulationMesh *parseCirclePacking(const std::string &vertexFile, const std::string &faceFile, const std::string &conformalFactorFile, double scale, double thickness);
    
    const MatrixX3m &getOriginalV() { return originalV; }
    const Eigen::MatrixX3i &getOriginalF() { return F; }

    void faceEnergyDensities(const MatrixX3m &V, VectorXm &densities);

    void testElasticEnergy();
    void lineSearch(const MatrixX3m &V, const std::vector<Matrix2m> &targetas, const VectorXm &deltaV, scalar &alpha, scalar &energy);
    void minimizeElasticEnergyNewton(MatrixX3m &V, int substeps, RenderingMesh &rm);
	void minimizeElasticEnergyNewtonWithMask(const MatrixX3m& mask, MatrixX3m &V, int substeps, RenderingMesh &rm);
	void saveObj(const MatrixX3m& V, const Eigen::MatrixX3i& F);
    std::vector<Matrix2m> baras;
	
	double shift;
private:
    void buildFaceWings();

	void randomShuffleDirection(VectorXm& deltaV, scalar shuffleRatio, scalar momentum);
	VectorXm maskFilter(const MatrixX3m & mask, const VectorXm & deltaV);

    MatrixX3m originalV;
	
    Eigen::MatrixX3i F;
    Eigen::MatrixX3i faceWings;

    VectorXm faceThicknesses;
	VectorXm randomVector;
	
    MaterialParameters params;
};

SimulationMesh *parseCirclePacking(const std::string &vertexFile, const std::string &faceFile, const std::string &conformalFactorFile, double scale, double thickness);

#endif
