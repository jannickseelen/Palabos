//
//  viscosity.cpp
//  Palabos
//
//  Created by MAC on 10/02/16.
//  Copyright © 2016 TUDelft. All rights reserved.
//

#include "palabos3D.h"
#include "palabos3D.hh"
#include <execinfo.h> // For complete stacktrace on exception see GNU libc manual
#include <iostream>

using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;

#define DESCRIPTOR descriptors::D3Q27Descriptor

class Var{
public:
	//Constructors
	Var(){
		plint extraLayer      = 0;              // Make the bounding box larger; for visualization purposes only. For the simulation, it is OK to have extraLayer=0.
		const plint blockSize = 0;             // Zero means: no sparse representation.
		const plint envelopeWidth = 1;          // For standard BGK dynamics.
		const plint extendedEnvelopeWidth = 2;  // Because the Guo off lattice boundary condition needs 2-cell neighbor access.
		bool performOutput = false;
		bool doImages = false;
		bool useAllDirections = false;
		bool useRegularizedWall = false;
		bool useIncompressible = false;
		bool poiseuilleInlet = false;
		bool convectiveScaling = false;

		T kinematicViscosity       = 0.;
		T averageInletVelocity     = 0.;
		plint referenceResolution  = 0;
		T nuLB                     = 0.;
		T fluidDensity             = 0.;
		T volume                   = 0.;
		T userDefinedInletDiameter = 0.;

		plint  = 0;
		plint openingSortDirection = 0;

		T simTime = 0;
		plint startLevel = 0;
		plint maxLevel   = 0;
		T epsilon = 0.;

		TriangleSet<T>* triangleSet = 0;
		TriangleSet<T>* triangleSetObstacle = 0;
		Array<T,3>* mount = 0;
		T currentTime = 0;
	}
	// Properties
	plint extraLayer, referenceResolution, , openingSortDirection, startLevel, maxLevel;
	const plint blockSize, envelopeWidth, extendedEnvelopeWidth;
	bool performOutput, doImages, useAllDirections, useRegularizedWall, useIncompressible, poiseuilleInlet, convectiveScaling;
	T kinematicViscosity, averageInletVelocity, nuLB, fluidDensity, volume, userDefinedInletDiameter, simTime, epsilon, currentTime;
	TriangleSet<T>* triangleSet, triangleSetObstacle;
	Array<T,3>* mount;
}
// Define Variables //


/* Obtain a backtrace and print it to stdout. */
void printTrace (void){
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);

  printf ("Obtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
     printf ("%s\n", strings[i]);

  free (strings);
}

void iniLattice( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,VoxelizedDomain3D<T>& voxelizedDomain ){
    // Switch all remaining outer cells to no-dynamics, except the outer
    //   boundary layer, and keep the rest as BGKdynamics.
    defineDynamics(lattice, voxelizedDomain.getVoxelMatrix(), lattice.getBoundingBox(),new NoDynamics<T,DESCRIPTOR>, voxelFlag::outside);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T) 1., Array<T,3>((T) 0.,(T) 0.,(T) 0.));
    lattice.initialize();
}

// This function outputs velocity, vorticity and pressure data, at selected
//   points of the computational domain, given their coordinates in physical units.
std::vector<T> pointMeasures (MultiBlockLattice3D<T,DESCRIPTOR>& lattice,Array<T,3> location, T dx, T dt ){
    std::vector<Array<T,3> > physicalPositions, positions;
    physicalPositions.push_back(Array<T,3>(0.022046, 0.015072, 0.044152));
    physicalPositions.push_back(Array<T,3>(0.027132, 0.049947, 0.095012));
    physicalPositions.push_back(Array<T,3>(0.034398, 0.056487, 0.057957));
    physicalPositions.push_back(Array<T,3>(0.031492, 0.025971, 0.084113));
    physicalPositions.push_back(Array<T,3>(0.025679, 0.025971, 0.091379));
    physicalPositions.push_back(Array<T,3>(0.018413, 0.011439, 0.076848));
    positions.resize(physicalPositions.size());

    for (pluint i=0; i<physicalPositions.size(); ++i) {
        positions[i] = (physicalPositions[i]-location)/dx;
    }

    std::vector<Array<T,3> > velocities = velocitySingleProbes(lattice, positions);
    std::vector<Array<T,3> > vorticities = vorticitySingleProbes(lattice, positions);
    std::vector<T> densities = densitySingleProbes(lattice, positions);

    std::vector<T> data;
    for (pluint i=0; i<physicalPositions.size(); ++i) {
        Array<T,3> pos = physicalPositions[i];
        Array<T,3> vel = velocities[i]*dx/dt;
        Array<T,3> vort = vorticities[i]/dt;
        T pressure = DESCRIPTOR<T>::cs2*(densities[i]-1.)*dx*dx/(dt*dt)*fluidDensity;
        if (performOutput) {
            pcout << "Pos ("
                  << pos[0] << "," << pos[1] << "," << pos[2]
                  << "); Velocity ("
                  << vel[0] << "," << vel[1] << "," << vel[2]
                  << "); Vorticity ("
                  << vort[0] << "," << vort[1] << "," << vort[2]
                  << "); Pressure " << pressure << std::endl;
        }
        data.push_back(norm(vel));
        data.push_back(norm(vort));
        data.push_back(pressure);
    }
    return data;
}


// Function to write data to image
void writeImages(OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>& boundaryCondition,Box3D const& imageDomain, 
Box3D const& vtkDomain, std::string fname, Array<T,3> location, T dx, T dt ){
    VtkImageOutput3D<T> vtkOut(fname, dx, location);
    vtkOut.writeData<float>(*boundaryCondition.computePressure(vtkDomain), "p", util::sqr(dx/dt)*fluidDensity);
    vtkOut.writeData<float>(*boundaryCondition.computeVelocityNorm(vtkDomain), "u", dx/dt);
    vtkOut.writeData<float>(*copyConvert<int,T>(*extractSubDomain(boundaryCondition.getVoxelizedDomain().getVoxelMatrix(),vtkDomain)), "voxel", 1.);

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledPpm(fname, *boundaryCondition.computeVelocityNorm(imageDomain));
}

// This function produces images at predefined yz, xz and xy planes. The coordinates of the planes are given
//   in physical coordinates, and the output variables are velocity, vorticity and pressure.
void writeImages(OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>& boundaryCondition, plint level, Array<T,3> location, T dx, T dt ){
    plint nx = boundaryCondition.getLattice().getNx();
    plint ny = boundaryCondition.getLattice().getNy();
    plint nz = boundaryCondition.getLattice().getNz();
    Array<T,3> yz_plane(0.016960, 0.032604, 0.057772);
    Array<T,3> xz_plane(0.026725, 0.017978, 0.057772);
    Array<T,3> xy_plane(0.026725, 0.032604, 0.084113);

    Array<T,3> lyz_plane((yz_plane-location)/dx);
    Array<T,3> lxz_plane((xz_plane-location)/dx);
    Array<T,3> lxy_plane((xy_plane-location)/dx);

    Box3D yz_imageDomain (
            util::roundToInt(lyz_plane[0]), util::roundToInt(lyz_plane[0]),
            0, ny-1, 0, nz-1 );
    Box3D xz_imageDomain (
            0, nx-1,
            util::roundToInt(lxz_plane[1]), util::roundToInt(lxz_plane[1]),
            0, nz-1 );
    Box3D xy_imageDomain (
            0, nx-1, 0, ny-1,
            util::roundToInt(lxy_plane[2]), util::roundToInt(lxy_plane[2]) );

    Box3D yz_vtkDomain (
            util::roundToInt(lyz_plane[0])-3, util::roundToInt(lyz_plane[0])+3,
            0, ny-1, 0, nz-1 );
    Box3D xz_vtkDomain (
            0, nx-1,
            util::roundToInt(lxz_plane[1])-3, util::roundToInt(lxz_plane[1])+3,
            0, nz-1 );
    Box3D xy_vtkDomain (
            0, nx-1, 0, ny-1,
            util::roundToInt(lxy_plane[2])-3, util::roundToInt(lxy_plane[2])+3 );

    writeImages(boundaryCondition, xy_imageDomain, xy_vtkDomain, "xy_"+util::val2str(level), location, dx, dt);
    writeImages(boundaryCondition, xz_imageDomain, xz_vtkDomain, "xz_"+util::val2str(level), location, dx, dt);
    writeImages(boundaryCondition, yz_imageDomain, yz_vtkDomain, "yz_"+util::val2str(level), location, dx, dt);
}

VoxelizedDomain3D<T> voxelizeBoundary(const bool& flowInside, TriangleBoundary3D<T>& boundary, const plint& borderWidth){
    if(flowInside){
		const int flowType = voxelFlag::inside; // For the outer boundaries
		const int borderType = voxelFlag::innerBorder; 
		} 
	else{
		const int flowType = voxelFlag::outside; // For the obstacle in the flow
		const int borderType = voxelFlag::outerBorder;
	} 
    VoxelizedDomain3D<T> voxelizedDomain(boundary, flowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);

	if (performOutput) {
        pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    }

    MultiScalarField3D<int> flagMatrix((MultiBlock3D&)voxelizedDomain.getVoxelMatrix());
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), flowType, flagMatrix.getBoundingBox(), 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(), borderType, flagMatrix.getBoundingBox(), 1);
    pcout << "Number of fluid cells: " << computeSum(flagMatrix) << std::endl;
	return voxelizedDomain;
}

	std::cout << "Creating Obstacle Structure" << '\n';
	// Move the obstacle to the starting point
    triangleSetObstacle.translate(mount);
	// Scale the triangle set
	DEFscaledMesh<T>* defMesh = new DEFscaledMesh<T>(*triangleSetObstacle, resolution, , margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
	// The immersed boundary method needs a set of vertices and a set of areas
    // that correspond to each vertex. These vertices and areas describe the
    // time dependent geometry of the surface at each time step.
    std::vector<Array<T,3> > vertices;
    std::vector<T> areas;
	// Get verticis and areas needed for immersed wall
    for (pluint iVertex = 0; iVertex < (pluint) defMesh->getMesh().getNumVertices(); iVertex++) {
        vertices.push_back(rectangleDef->getMesh().getVertex(iVertex));
        areas.push_back(rectangleDef->getMesh().computeVertexArea(iVertex));
    }
	delete defMesh;
	
    boundary.getMesh().inflate();
	// When convective scaling is used (relationship of dt with respect to dx as the grid is
    //   refined) the value of the kinematic viscosity must be also properly adjusted.
    T nuLB_ = nuLB;
    if (convectiveScaling) {
        nuLB_ = nuLB * util::twoToThePower(level);
    }
    T dx = boundary.getDx();
    T dt = nuLB_ / kinematicViscosity *dx*dx;
    T omega = 1./(3.*nuLB_+0.5);
    Array<T,3> location(boundary.getPhysicalLocation());

    if (performOutput) {
		pcout << "nuLB=" << nuLB_ << std::endl;
		pcout << "tau=" << 1./omega << std::endl;
        pcout << "dx=" << dx << std::endl;
        pcout << "dt=" << dt << std::endl;
    }
	VoxelizedDomain3D<T> voxelizedDomain = voxelizeBoundary(false, boundary, borderWidth);
}

// This is the function that prepares and performs the actual simulation. (NOT FINISHED YET)
std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > run (plint level, MultiBlockLattice3D<T,DESCRIPTOR>* iniVal=0 ){
    plint margin = 3;           // Extra margin of allocated cells around the obstacle.
    plint borderWidth = 1;      // Because the Guo boundary condition acts in a one-cell layer.
    if(margin<borderWidth)
		throw "Margin cannot be smaller then borderWidth";	// Requirement: margin>=borderWidth.

    // The resolution is doubled at each coordinate direction with the increase of the
    //   resolution level by one. The parameter ``referenceResolution'' is by definition
    //   the resolution at grid refinement level 0.
    plint resolution = referenceResolution * util::twoToThePower(level);

    // The next few lines of code are typical. They transform the surface geometry of the
    //   mesh given by the user to more efficient data structures that are internally
    //   used by palabos. The TriangleBoundary3D structure will be later used to assign
    //   proper boundary conditions.
    DEFscaledMesh<T>* defMesh = new DEFscaledMesh<T>(*triangleSet, resolution, , margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();

    // When convective scaling is used (relationship of dt with respect to dx as the grid is
    //   refined) the value of the kinematic viscosity must be also properly adjusted.
    T nuLB_ = nuLB;
    if (convectiveScaling) {
        nuLB_ = nuLB * util::twoToThePower(level);
    }
    T dx = boundary.getDx();
    T dt = nuLB_ / kinematicViscosity *dx*dx;
    T omega = 1./(3.*nuLB_+0.5);
    Array<T,3> location(boundary.getPhysicalLocation());

    pcout << "nuLB=" << nuLB_ << std::endl;
    pcout << "tau=" << 1./omega << std::endl;
    if (performOutput) {
        pcout << "dx=" << dx << std::endl;
        pcout << "dt=" << dt << std::endl;
    }

	VoxelizedDomain3D<T> voxelizedDomain = voxelizeBoundary(true, boundary, borderWidth);

    Dynamics<T,DESCRIPTOR>* dynamics = 0;
    if (useIncompressible) { dynamics = new IncBGKdynamics<T,DESCRIPTOR>(omega); }// In this model velocity equals momentum.
    else {dynamics = new BGKdynamics<T,DESCRIPTOR>(omega); }// In this model velocity equals momentum divided by density.
	
    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice 
	= generateMultiBlockLattice<T,DESCRIPTOR> (voxelizedDomain.getVoxelMatrix(), envelopeWidth, dynamics );
    lattice->toggleInternalStatistics(false);

    // The next piece of code is put for efficiency reasons at communications in parallel runs.
    //   The efficiency advantage comes essentially because the density and velocity are
    //   written in different fields.
    std::vector<MultiBlock3D*> rhoBarJarg;
    plint numScalars = 4;
    MultiNTensorField3D<T>* rhoBarJfield =
    generateMultiNTensorField3D<T>(*lattice, extendedEnvelopeWidth, numScalars);
    rhoBarJfield->toggleInternalStatistics(false);
    rhoBarJarg.push_back(rhoBarJfield);
    plint processorLevel=0;
    integrateProcessingFunctional (new PackedRhoBarJfunctional3D<T,DESCRIPTOR>(),lattice->getBoundingBox(), *lattice, *rhoBarJfield, processorLevel );

    // The Guo off lattice boundary condition is set up.
	const int flowType = voxelFlag::inside;
	BoundaryProfiles3D<T,Velocity> profiles;
    OffLatticeModel3D<T,DESCRIPTOR>* model =
    new OffLatticeModel3D<T,DESCRIPTOR> (
		new TriangleFlowShape3D<T,Array<T,3> > (voxelizedDomain.getBoundary(), profiles),flowType, useAllDirections );
    model->setVelIsJ(useIncompressible); // When the incompressible BGK model is used, velocity equals momentum.
    model->selectUseRegularizedModel(useRegularizedWall);
    model->selectComputeStat(false);
    OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity> boundaryCondition (model, voxelizedDomain, *lattice);
    boundaryCondition.insert(rhoBarJarg);

    iniLattice(*lattice, voxelizedDomain);
    if(iniVal) {
        Box3D toDomain(lattice->getBoundingBox());
        Box3D fromDomain(toDomain.shift(margin,margin,margin)); // During rescaling, the margin doubled in size,
        //   an effect which is cancelled here through a shift.
        copy(*iniVal, fromDomain, *lattice, toDomain, modif::staticVariables);
        computePackedRhoBarJ(*lattice, *rhoBarJfield, lattice->getBoundingBox());
        boundaryCondition.apply(rhoBarJarg);
    }

    // The ValueTracer is needed to check when a chosen quantity (in our case the average energy)
    //   has converged, so to conclude that steady state has been reached for the specific grid
    //   refinement level and stop the simulation.
    plint convergenceIter=20;
    util::ValueTracer<T> velocityTracer(0.05*convergenceIter, resolution, epsilon);
    global::timer("iteration").restart();
    plint i = util::roundToInt(currentTime/dt);
    lattice->resetTime(i);

	// The next container block is necessary for the immersed-wall algorithm.
    MultiContainerBlock3D container(*rhoBar);

    // Collision and streaming iterations.
    while(!velocityTracer.hasConverged() && currentTime<simTime)
    {
        if (i%200==0 && performOutput) {
            pcout << "T= " << currentTime << "; "
            << "Average energy: "
            << boundaryCondition.computeAverageEnergy()*util::sqr(dx/dt) << std::endl;
        }
        if (i%convergenceIter==0) {
            velocityTracer.takeValue(computeAverageEnergy(*lattice));
        }
		
		lattice->executeInternalProcessors(); // Execute all processors and communicate appropriately.
		
		// Immersed walls algorithm.
        T timeLB = currentTime + 1.0;
		
		// Instantiate the immersed wall data and performed the immersed boundary iterations.
        instantiateImmersedWallData(vertices, areas, container);
		
		// see offLattice / immersedWall3D.hh
		computeImmersedBoundaryForce();
		
		
        lattice->collideAndStream();

        ++i;
        currentTime = i*dt;
    }
    delete rhoBarJfield;

    Box3D measureBox(lattice->getBoundingBox());

    // Image output.
    if (doImages) {
        writeImages(boundaryCondition, level, location, dx, dt);
        std::vector<std::string> scalarNames;
        scalarNames.push_back("pressure");
        scalarNames.push_back("wss");
        std::vector<T> scalarFactor;
        scalarFactor.push_back(util::sqr(dx/dt)*fluidDensity);
        scalarFactor.push_back(util::sqr(dx/dt)*fluidDensity);
        std::vector<std::string> vectorNames;
        vectorNames.push_back("force");
        std::vector<T> vectorFactor;
        vectorFactor.push_back(util::sqr(dx/dt)*fluidDensity);
        bool dynamicMesh = false;
        writeSurfaceVTK (
                         boundary,
                         *computeSurfaceForce( boundary, voxelizedDomain, *lattice, model->velIsJ(), dynamicMesh ),
                         scalarNames, vectorNames, "surface_"+util::val2str(level)+".vtk", dynamicMesh, 0,
                         scalarFactor, vectorFactor );
    }

    T averageEnergy = boundaryCondition.computeAverageEnergy()*util::sqr(dx/dt);
    T rmsVorticity  = boundaryCondition.computeRMSvorticity()/dt;
    T pressureDrop = util::sqr(dx/dt)*fluidDensity;

    if (performOutput) {
        pcout << "Average energy: " << averageEnergy << std::endl;
        pcout << "Total energy: " << averageEnergy*volume << std::endl;
        pcout << "RMS vorticity * volume * 0.5: " << rmsVorticity*0.5*volume << std::endl;
        pcout << "Pressure drop: " << pressureDrop << std::endl;
        pcout << "Number of iterations: " << i << std::endl;
    }
    pcout << "Elapsed time: " << global::timer("iteration").stop() << std::endl;
    pcout << "Total elapsed time: " << global::timer("global").getTime() << std::endl;

    if (performOutput) {
        pcout << "Description: "
        << "Tot. energy, pressure-drop, tot. enstrophy,"
        << "  vel1, vort1, pres1,  vel2, vort2, pres2,"
        << "  vel3, vort3, pres3,  vel4, vort4, pres4,"
        << "  vel5, vort5, pres5,  vel6, vort6, pres6" << std::endl;
        pcout << "All data: ";
    }
    pcout << averageEnergy*volume << ", " << pressureDrop << ", " << rmsVorticity*volume*0.5 << ", ";
    std::vector<T> pointData = pointMeasures(*lattice, location, dx, dt);
    for (pluint i=0; i<pointData.size(); ++i) {
        pcout << pointData[i];
        if (i!=pointData.size()-1) {
            pcout << ", ";
        }
    }
    pcout << std::endl;
    
    return lattice;
}

// Read the user input XML file//
void readParameters(XMLreader const& document){
    std::string meshFileName;
	std::string meshFileNameObstacle;
    document["geometry"]["mesh"].read(meshFileName);
	document["geometry"]["obstacle"].read(meshFileNameObstacle);
	
	double mountX; double mountY; double mountZ;
	document["mount"]["x"].read(mountX);
	document["mount"]["y"].read(mountY);
	document["mount"]["z"].read(mountZ); 

    document["fluid"]["kinematicViscosity"].read(kinematicViscosity);
    document["fluid"]["density"].read(fluidDensity);
    document["fluid"]["volume"].read(volume);

    document["numerics"][""].read();
    document["numerics"]["referenceResolution"].read(referenceResolution);
    document["numerics"]["nuLB"].read(nuLB);

    document["simulation"]["simTime"].read(simTime);
    document["simulation"]["maxLevel"].read(maxLevel);
    document["simulation"]["epsilon"].read(epsilon);

    document["simulation"]["performOutput"].read(performOutput);
    document["simulation"]["doImages"].read(doImages);
    document["simulation"]["useAllDirections"].read(useAllDirections);
    document["simulation"]["useRegularizedWall"].read(useRegularizedWall);
    document["simulation"]["useIncompressible"].read(useIncompressible);
    document["simulation"]["poiseuilleInlet"].read(poiseuilleInlet);
    document["simulation"]["convectiveScaling"].read(convectiveScaling);

    // At this part, the surface geometry (as given by the user in
    //   the form of an ASCII or binary STL file) is read into a data structure
    //   comprised by a set of triangles. The DBL constant means that double
    //   precision accuracy will be used (generally the recommended choice).
    triangleSet = new TriangleSet<T>(meshFileName, DBL);
	traingleSetObstacle = new TriangleSet<T>(meshFileNameObstacle, DBL);
	
	// Set mount point of the obstacle
	mount = new Array<T,3>(mountX,mountY,mountZ);
}
	
int xmlProcess(char* argv[]){
	string paramXmlFileName;
    try {global::argv(1).read(paramXmlFileName);}
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: "
        << (std::string)global::argv(0) << " parameter-input-file.xml" << std::endl;
        return -1;
    }
    // Read the parameter XML input file. (Lots of comments are included there too).
    try {
        XMLreader document(paramXmlFileName);
        readParameters(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Error in input file " << paramXmlFileName
        << ": " << exception.what() << std::endl;
        return -1;
    }
	if(referenceResolution==0){
		throw "ERROR: Resolution cannot be zero!";
	}
	return 0;
}

int main(int argc, char* argv[]) {
	plb::plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    global::IOpolicy().activateParallelIO(false);

    if(xmlProcess(argv) != 0){return -1;}
	
    global::timer("global").start();
    plint iniLevel=0;
    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > iniConditionLattice(0);
    // This code incorporates the concept of smooth grid refinement until convergence is
    //   achieved. The word ``smooth'' indicates that as the refinement level increases
    //   by one, the whole grid doubles in each direction. When the grid is refined, both
    //   dx and dt have to change. Whether dt is changed as dx^2 (diffusive behavior)
    //   or as dx (convective behavior), is controlled by the input variable
    //   ``convectiveScaling'' (the recommended choice is not to use convective scaling).
    try {
        for (plint level=iniLevel; level<=maxLevel; ++level) {
            pcout << std::endl << "Running new simulation at level " << level << std::endl;
            std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > convergedLattice (run(level, iniConditionLattice.get()) );
            if (level != maxLevel) {
                plint dxScale = -1;
                plint dtScale = -2;
                if (convectiveScaling) {dtScale = -1;}
                // The converged simulation of the previous grid level is used as the initial condition
                //   for the simulation at the next grid level (after appropriate interpolation has
                //   taken place).
                iniConditionLattice = std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > (refine(*convergedLattice, dxScale, dtScale, new BGKdynamics<T,DESCRIPTOR>(1.)) );
            }
        }
    }
    catch(PlbException& exception) {
		printTrace();
        pcout << exception.what() << std::endl;
        return -1;
    }
    return 0;
}
