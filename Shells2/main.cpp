#ifdef GUI_VERSION
#include <igl/viewer/Viewer.h>
#endif
#include "SimulationMesh.h"
#include "RenderingMesh.h"
#ifdef GUI_VERSION
#include <thread>
#endif
#include <igl/writeOBJ.h>
#include "ArbitraryPrecision.h"
#include <fstream>
#include <sstream>
#include <string>

SimulationMesh *sm;
MaterialParameters params;
RenderingMesh rm;
scalar shift0;
scalar shift1;
scalar shiftstep = 0.00001;
int processNo;

void perturbStartV(MatrixX3m & startV, MatrixX3m & mask, scalar shift)
{
	int nverts = startV.rows();
	mask.resize(nverts, 3);
	mask.setOnes();

	int edgeLength = rint(sqrt(nverts));

	for (int i = 0; i < edgeLength; ++i)
	{
		int index0 = i * edgeLength;
		int index1 = index0 + edgeLength - 1;

	//	startV.row(index0) += Vector3m(0, shift, 0);
	//	startV.row(index1) += Vector3m(0, -shift, 0);

		mask.row(index0).setZero();
		mask.row(index1).setZero();
	}

	for (int i = 0; i < nverts; ++i)
	{
		startV(i, 1) = 0.5 + (1 - 2 * shift) * (startV(i, 1) - 0.5);
	}

}

void saveV(std::ofstream& fout, const MatrixX3m& V)
{
	for (int i = 0; i < V.rows(); ++i)
	{
		fout << V(i, 0) << ' '
		<< V(i, 1) << ' '
		<< V(i, 2) << ' ';
	}
	fout << std::endl;
}

void optimize()
{
	std::stringstream ss;
	ss << "square_8x8_shrink_init_" << processNo << ".txt";
		  
	std::ofstream initOut(ss.str(), std::ios::app);
	
	std::stringstream ss1;
	ss1 << "square_8x8_shrink_real_" << processNo << ".txt";

	std::ofstream realOut(ss1.str(), std::ios::app);
	//std::ofstream maskOut("square_8x8_shrink_mask_" + str + ".txt", std::ios::app);

	for (scalar s = shift0; s < shift1; s = s + shiftstep)
	{
    		MatrixX3m startV = sm->getOriginalV();
		MatrixX3m mask;
		perturbStartV(startV, mask, s);
		MatrixX3m startV0(startV);
		//saveV(initOut, startV);
		//saveV(maskOut, mask);

   		//sm->minimizeElasticEnergyNewton(startV, 20, rm);
		sm->minimizeElasticEnergyNewtonWithMask(mask, startV, 1, rm);
		saveV(initOut, startV0);
		//saveV(maskOut, mask);
		saveV(realOut, startV);
	}
	std::cout << "Optimization done" << std::endl;
}

#ifdef GUI_VERSION
bool pre_draw(igl::viewer::Viewer &viewer)
{
    viewer.data.clear();
    rm.visMutex.lock();
    rm.recomputeColors();
    viewer.data.set_mesh(rm.V, rm.F);
    viewer.data.set_colors(rm.C);
    
    int ngrads = rm.gradientDir.rows();
    Eigen::MatrixXd P1(ngrads/3, 3);
    Eigen::MatrixXd P2(ngrads/3, 3);
    Eigen::MatrixXd eC(ngrads/3, 3);
    eC.setOnes();
    for (int i = 0; i < ngrads; i++)
    {
        P1(i / 3, i % 3) = rm.V(i / 3, i % 3);
        P2(i / 3, i % 3) = rm.V(i / 3, i % 3) - rm.gradScale*rm.gradientDir[i];
        
    }

    int ndescent = rm.descentDir.rows();
    Eigen::MatrixXd D1(ndescent / 3, 3);
    Eigen::MatrixXd D2(ndescent / 3, 3);
    Eigen::MatrixXd dC(ndescent / 3, 3);
    dC.setZero();
    for (int i = 0; i < ndescent; i++)
    {        
        D1(i / 3, i % 3) = rm.V(i / 3, i % 3);
        D2(i / 3, i % 3) = rm.V(i / 3, i % 3) - rm.gradScale*rm.descentDir[i];

    }
    rm.visMutex.unlock();

    if(rm.showGradientDir)
        viewer.data.add_edges(P1, P2, eC);
    if(rm.showDescentDir)
        viewer.data.add_edges(D1, D2, dC);
    viewer.data.set_face_based(true);
    return false;
}

// This function is called every time a keyboard button is pressed
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
    if (key == 'g' || key == 'G')
    {
        rm.showGradientDir = !rm.showGradientDir;
    }
    if (key == 's' || key == 'S')
    {
        std::thread t(optimize);
        t.detach();
    }
    else if (key == '[')
    {
        rm.colorMax /= 2.0;
        std::cout << "New max color: " << rm.colorMax << std::endl;        
    }
    else if (key == ']')
    {
        rm.colorMax *= 2.0;
        std::cout << "New max color: " << rm.colorMax << std::endl;
    }
    else if (key == ',')
    {
        rm.gradScale /= 2.0;
        std::cout << "New grad scale: " << rm.gradScale << std::endl;
    }
    else if (key == '.')
    {
        rm.gradScale *= 2.0;
        std::cout << "New grad scale: " << rm.gradScale << std::endl;
    }
    else if (key == 'd' || key == 'D')
    {
        rm.showDescentDir = !rm.showDescentDir;
    }
    return true;
}
#endif

// Rountines for constrained cloth
int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		return -1;
	}
	
	//int label0 = atoi(argv[1]);
	//int label1 = atoi(argv[2]);
	processNo = atoi(argv[1]);

	//shift0 = label0 * 0.001;
	//shift1 = label1 * 0.001;
	shift0 = processNo * 0.001;
	shift1 = (processNo + 1) * 0.001;

	//mpfr::mpreal::set_default_prec(256);

	std::ifstream streamfile("0_square.txt");
	if (!streamfile)
	{
		std::cerr << "Couldn't open settings file " << argv[1] << std::endl;
		return -1;
	}

	std::string verts, faces, conffactors, outpath;
	getline(streamfile, verts);
	getline(streamfile, faces);
	getline(streamfile, conffactors);
	getline(streamfile, outpath);

	sm = parseCirclePacking("centers_0_square.dat", "faces_0_square.dat", "conformal_0_square.dat", 1.0, 0.05);
	//sm->shift = double(shift);
	/*Eigen::MatrixX3d V(4, 3);
	V << 0, 0, 0,
	1, 0, 0,
	0, 1, 0,
	-1, 0, 0;
	Eigen::MatrixX3i F(2, 3);
	F << 0, 1, 2,
	2, 3, 0;

	sm = new SimulationMesh(V, F, 1e-1);
	sm->baras[0] *= 2;
	sm->baras[1] *= 1.0;*/
#ifdef GUI_VERSION
	rm.outputPath = outpath;
	if (!sm)
	{
		std::cerr << "Couldn't load input files" << std::endl;
		return -1;
	}

	int nverts = (int)sm->getOriginalV().rows();

	rm.V.resize(sm->getOriginalV().rows(), 3);
	for (int i = 0; i<(int)sm->getOriginalV().rows(); i++)
		for (int j = 0; j<3; j++)
			rm.V(i, j) = double(sm->getOriginalV()(i, j));
	rm.F = sm->getOriginalF();
	rm.C.resize(rm.F.rows(), 3);
	rm.C.setZero();

	//sm->testElasticEnergy();    
	//std::cout << rm.V << std::endl;

	igl::viewer::Viewer viewer;
	viewer.core.is_animating = true;
	viewer.callback_key_down = &key_down;
	viewer.callback_pre_draw = &pre_draw;
	viewer.launch();
#else
	optimize();
#endif
}


// Original Main Function
//int main(int argc, char *argv[])
//{
//    if (argc != 2)
//        return -1;
//
//    mpfr::mpreal::set_default_prec(256);
//
//    std::ifstream streamfile(argv[1]);
//    if (!streamfile)
//    {
//        std::cerr << "Couldn't open settings file " << argv[1] << std::endl;
//        return -1;
//    }
//
//    std::string verts, faces, conffactors, outpath;
//    getline(streamfile, verts);
//    getline(streamfile, faces);
//    getline(streamfile, conffactors);
//    getline(streamfile, outpath);
//
//    sm = parseCirclePacking(verts, faces, conffactors, 1.0, 1e-1);
//    /*Eigen::MatrixX3d V(4, 3);
//    V << 0, 0, 0,
//        1, 0, 0,
//        0, 1, 0,
//        -1, 0, 0;
//    Eigen::MatrixX3i F(2, 3);
//    F << 0, 1, 2,
//        2, 3, 0;
//
//    sm = new SimulationMesh(V, F, 1e-1);
//    sm->baras[0] *= 2;
//    sm->baras[1] *= 1.0;*/
//    rm.outputPath = outpath;    
//    if (!sm)
//    {
//        std::cerr << "Couldn't load input files" << std::endl;
//        return -1;
//    }
//
//    int nverts = (int)sm->getOriginalV().rows();
//    
//    rm.V.resize(sm->getOriginalV().rows(), 3);
//    for(int i=0; i<(int)sm->getOriginalV().rows(); i++)
//        for(int j=0; j<3; j++)
//            rm.V(i,j) = double(sm->getOriginalV()(i,j));
//    rm.F = sm->getOriginalF();
//    rm.C.resize(rm.F.rows(), 3);
//    rm.C.setZero();    
//
//    //sm->testElasticEnergy();    
//    //std::cout << rm.V << std::endl;
//
//    igl::viewer::Viewer viewer;
//    viewer.core.is_animating = true;
//    viewer.callback_key_down = &key_down;
//    viewer.callback_pre_draw = &pre_draw;
//    viewer.launch();   
//}
