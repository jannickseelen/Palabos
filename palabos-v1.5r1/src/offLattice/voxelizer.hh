/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: Orestis Malaspinas */

#ifndef VOXELIZER_HH
#define VOXELIZER_HH

#include "core/globalDefs.h"
#include "offLattice/voxelizer.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include <dataProcessors/metaStuffWrapper3D.hh>
#include <dataProcessors/dataInitializerWrapper3D.hh>
#include <core/plbTimer.h>
#include <thread>
#include <modules/mpiDataManager.hh>
#include <memory>
#include <vector>
#include <cmath>

namespace plb {

namespace voxelFlag {
    inline int invert(int arg) {
        switch(arg) {
            case inside: return outside;
            case outside: return inside;
            case innerBorder: return outerBorder;
            case outerBorder: return innerBorder;
            case undetermined: return undetermined;
            default:
                PLB_ASSERT(false);
        }
        return undetermined;
    }
    inline int bulkFlag(int arg) {
        if (arg==innerBorder || arg==inside) {
            return inside;
        }
        else if (arg==outerBorder || arg==outside) {
            return outside;
        }
        else {
            return undetermined;
        }
    }
    inline int borderFlag(int arg) {
        if (arg==inside || arg==innerBorder) {
            return innerBorder;
        }
        else if (arg==outside || arg==outerBorder) {
            return outerBorder;
        }
        else {
            return undetermined;
        }
    }
    inline bool insideFlag(int arg) {
        return arg==inside || arg==innerBorder;
    }
    inline bool outsideFlag(int arg) {
        return arg==outside || arg==outerBorder;
    }

}  // namespace voxelFlag

template<typename T>
std::unique_ptr<MultiScalarField3D<int> > voxelize (
        TriangularSurfaceMesh<T> const& mesh,
        plint symmetricLayer, plint borderWidth )
{
    Array<T,2> xRange, yRange, zRange;
    mesh.computeBoundingBox(xRange, yRange, zRange);
    // Creation of the multi-scalar field. The +1 is because if the resolution is N,
    //   the number of nodes is N+1.
    plint nx = (plint)(xRange[1] - xRange[0]) + 1 + 2*symmetricLayer;
    plint ny = (plint)(yRange[1] - yRange[0]) + 1 + 2*symmetricLayer;
    plint nz = (plint)(zRange[1] - zRange[0]) + 1 + 2*symmetricLayer;
	std::unique_ptr<MultiScalarField3D<int> > pointer = voxelize(mesh, Box3D(0,nx-1, 0,ny-1, 0,nz-1), borderWidth);
    return pointer;
}

template<typename T>
std::unique_ptr<MultiScalarField3D<int> > voxelize (TriangularSurfaceMesh<T> const& mesh, Box3D const& domain, plint borderWidth )
{
	#ifdef PLB_DEBUG
		bool main = false;
		main = global::mpi().isMainProcessor();
		if(main){std::cout<< "[DEBUG] Voxelizing Triangular Surface Mesh"<<std::endl;}
	#endif
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    plint envelopeWidth=1;
    std::unique_ptr<MultiScalarField3D<int> > voxelMatrix = generateMultiScalarField<int>(domain, voxelFlag::outside, envelopeWidth);
	Box3D voxelDomain = voxelMatrix->getBoundingBox();

    setToConstant(*voxelMatrix, voxelDomain.enlarge(-1), voxelFlag::undetermined);
    MultiContainerBlock3D hashContainer(*voxelMatrix);
    std::vector<MultiBlock3D*> container_arg;
    container_arg.push_back(&hashContainer);

    applyProcessingFunctional(new CreateTriangleHash<T>(mesh),hashContainer.getBoundingBox(), container_arg);

    std::vector<MultiBlock3D*> flag_hash_arg;
    flag_hash_arg.push_back(voxelMatrix.get());
    flag_hash_arg.push_back(&hashContainer);

    voxelMatrix->resetFlags(); // Flags are used internally by VoxelizeMeshFunctional3D.
	while (!allFlagsTrue(voxelMatrix.get())) {
		VoxelizeMeshFunctional3D<T>* func = new VoxelizeMeshFunctional3D<T>(mesh);
		applyProcessingFunctional(func, voxelDomain, flag_hash_arg);
	}

	detectBorderLine(*voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth);
	#ifdef PLB_DEBUG
		if(main){std::cout<< "[DEBUG] DONE Voxelizing Triangular Surface Mesh"<<std::endl;}
	#endif
    return voxelMatrix;
}

template<typename T>
std::unique_ptr<MultiScalarField3D<int> > voxelize (
        TriangularSurfaceMesh<T> const& mesh,
        Box3D const& domain, plint borderWidth, Box3D seed )
{
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    plint envelopeWidth=1;

    std::unique_ptr<MultiScalarField3D<int> > voxelMatrix
        = generateMultiScalarField<int>(domain, voxelFlag::undetermined, envelopeWidth);
    setToConstant(*voxelMatrix, seed, voxelFlag::outside);

    MultiContainerBlock3D hashContainer(*voxelMatrix);
    std::vector<MultiBlock3D*> container_arg;
    container_arg.push_back(&hashContainer);
    applyProcessingFunctional (
            new CreateTriangleHash<T>(mesh),
            hashContainer.getBoundingBox(), container_arg );

    std::vector<MultiBlock3D*> flag_hash_arg;
    flag_hash_arg.push_back(voxelMatrix.get());
    flag_hash_arg.push_back(&hashContainer);

    voxelMatrix->resetFlags(); // Flags are used internally by VoxelizeMeshFunctional3D.
    plint maxIteration=100;
    plint i=0;
    while (!allFlagsTrue(voxelMatrix.get()) && i<maxIteration) {
        applyProcessingFunctional (
                new VoxelizeMeshFunctional3D<T>(mesh),
                voxelMatrix->getBoundingBox(), flag_hash_arg );
        ++i;
    }
    if (i==maxIteration) {
        pcout << "Warning: Voxelization failed." << std::endl;
    }

    detectBorderLine(*voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth);

    return voxelMatrix;
}


template<typename T>
std::unique_ptr<MultiScalarField3D<int> > revoxelize (TriangularSurfaceMesh<T> const& mesh, MultiScalarField3D<int>& oldVoxelMatrix,
MultiContainerBlock3D& hashContainer, plint borderWidth ){
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    Box3D domain(oldVoxelMatrix.getBoundingBox());
    std::unique_ptr<MultiScalarField3D<int> > voxelMatrix (
            new MultiScalarField3D<int>((MultiBlock3D&)oldVoxelMatrix) );
    setToConstant(*voxelMatrix, domain, voxelFlag::outside);
    setToConstant(*voxelMatrix, voxelMatrix->getBoundingBox().enlarge(-1),
                  voxelFlag::undetermined);

    std::vector<MultiBlock3D*> flag_hash_arg;
    flag_hash_arg.push_back(voxelMatrix.get());
    flag_hash_arg.push_back(&hashContainer);

    voxelMatrix->resetFlags(); // Flags are used internally by VoxelizeMeshFunctional3D.
    while (!allFlagsTrue(voxelMatrix.get())) {
        applyProcessingFunctional (
                new VoxelizeMeshFunctional3D<T>(mesh),
                voxelMatrix->getBoundingBox(), flag_hash_arg );
    }

    detectBorderLine(*voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth);

    return voxelMatrix;
}


/* ******** VoxelizeMeshFunctional3D ************************************* */

template<typename T>
VoxelizeMeshFunctional3D<T>::VoxelizeMeshFunctional3D (TriangularSurfaceMesh<T> const& mesh_): mesh(mesh_){ }

template<typename T>
bool VoxelizeMeshFunctional3D<T>::distanceToSurface (AtomicContainerBlock3D& hashContainer, Array<T,3> const& point,
 T& distance, bool& isBehind) const{
    T maxDistance = std::sqrt((T)3);
    Array<T,2> xRange(point[0]-maxDistance, point[0]+maxDistance);
    Array<T,2> yRange(point[1]-maxDistance, point[1]+maxDistance);
    Array<T,2> zRange(point[2]-maxDistance, point[2]+maxDistance);
    TriangleHash<T> triangleHash(hashContainer);
    std::vector<plint> possibleTriangles;
    triangleHash.getTriangles(xRange, yRange, zRange, possibleTriangles);

    T    tmpDistance;
    bool tmpIsBehind;
    bool triangleFound = false;

    for (pluint iPossible=0; iPossible<possibleTriangles.size(); ++iPossible) {
        plint iTriangle = possibleTriangles[iPossible];
        mesh.distanceToTriangle (
                    point, iTriangle, tmpDistance, tmpIsBehind );
        if( !triangleFound || tmpDistance<distance) {
            distance = tmpDistance;
            isBehind = tmpIsBehind;
            triangleFound = true;
        }
    }
    return triangleFound;
}

template<typename T>
bool VoxelizeMeshFunctional3D<T>::checkIfFacetsCrossed (AtomicContainerBlock3D& hashContainer,
Array<T,3> const& point1, Array<T,3> const& point2, T& distance, plint& whichTriangle)
{
    Array<T,2> xRange (
                 std::min(point1[0], point2[0]),
                 std::max(point1[0], point2[0]) );
    Array<T,2> yRange (
                 std::min(point1[1], point2[1]),
                 std::max(point1[1], point2[1]) );
    Array<T,2> zRange (
                 std::min(point1[2], point2[2]),
                 std::max(point1[2], point2[2]) );
    TriangleHash<T> triangleHash(hashContainer);
    std::vector<plint> possibleTriangles;
    triangleHash.getTriangles(xRange, yRange, zRange, possibleTriangles);

    int flag = 0; // Check for crossings inside the point1-point2 segment.
    Array<T,3> intersection; // Dummy variable.
    Array<T,3> normal;       // Dummy variable.
    T tmpDistance;           // Dummy variable.

    if (global::counter("voxelizer-debug").getCount()==1) {
        std::cout << "{";
    }
    std::vector<T> crossings;
    for (pluint iPossible=0; iPossible<possibleTriangles.size(); ++iPossible) {
        plint iTriangle = possibleTriangles[iPossible];
        if (mesh.pointOnTriangle(point1, point2, flag, iTriangle, intersection, normal, tmpDistance)==1) {
            if (global::counter("voxelizer-debug").getCount()==1) {
                std::cout << "(" << iTriangle << ";" << tmpDistance << ")";
            }
            crossings.push_back(tmpDistance);
            if (crossings.size()==1 || tmpDistance<distance) {
                distance = tmpDistance;
                whichTriangle = iTriangle;
            }
        }
    }
    if (global::counter("voxelizer-debug").getCount()==1) {
        std::cout << "}";
    }

    if (crossings.size()==0) {
        return false;
    }
    else {
        bool hasCrossed = true;
        for (pluint iCrossing=1; iCrossing<crossings.size(); ++iCrossing) {
            //const T eps1 = std::numeric_limits<double>::epsilon()*1.e2;
            //if ( !util::fpequal(crossings[iCrossing], crossings[iCrossing-1], eps1) )

            //const T eps1 = std::numeric_limits<double>::epsilon()*1.e4;
            
            const T eps1 = std::numeric_limits<double>::epsilon()*1.e4;
            if ( std::fabs(crossings[iCrossing]-crossings[iCrossing-1])>eps1)
            {
                hasCrossed = !hasCrossed;
            }
        }
        return hasCrossed;
    }
}

template<typename T>
bool VoxelizeMeshFunctional3D<T>::createVoxelizationRange (Box3D const& domain, ScalarField3D<int>& voxels, Array<plint,2>& xRange,
Array<plint,2>& yRange, Array<plint,2>& zRange ){
    // The purpose of the first three loops is to locate the eight
    //   corners of the cube. One voxel per corner would be insufficient
    //   because a potential seed is situated differently, depending on
    //   whether it is on the boundary of the multi-block or somewhere inside.
    for (plint dx=0; dx<=+1; ++dx) {
        plint xMin = domain.x0+dx*domain.getNx()-1;
        plint xMax = domain.x0+dx*domain.getNx();
        for (plint dy=0; dy<=+1; ++dy) {
            plint yMin = domain.y0+dy*domain.getNy()-1;
            plint yMax = domain.y0+dy*domain.getNy();
            for (plint dz=0; dz<=+1; ++dz) {
                plint zMin = domain.z0+dz*domain.getNz()-1;
                plint zMax = domain.z0+dz*domain.getNz();

                // Locate a potential seed in one of the corners.
                for (plint iX=xMin; iX<=xMax; ++iX) {
                    for (plint iY=yMin; iY<=yMax; ++iY) {
                        for (plint iZ=zMin; iZ<=zMax; ++iZ) {
                            if (voxels.get(iX,iY,iZ) != voxelFlag::undetermined) {
                                xRange[0] = domain.x0+dx*(domain.getNx()-1);
                                xRange[1] = domain.x0+(1-dx)*(domain.getNx()-1);
                                yRange[0] = domain.y0+dy*(domain.getNy()-1);
                                yRange[1] = domain.y0+(1-dy)*(domain.getNy()-1);
                                zRange[0] = domain.z0+dz*(domain.getNz()-1);
                                zRange[1] = domain.z0+(1-dz)*(domain.getNz()-1);
                                return true;
                            }
                        }
                    }
                }

            }
        }
    }
    return false;
}

template<typename T>
void VoxelizeMeshFunctional3D<T>::printOffender (ScalarField3D<int> const& voxels, AtomicContainerBlock3D& hashContainer,
const Dot3D& pos){
    std::set<plint> triangles;
    Dot3D offset = voxels.getLocation();
    Dot3D pos_ = pos+offset;
    std::cout << "Position (" << pos_.x << "," << pos_.y << "," << pos_.z << ")" << std::endl;
    for (plint dx=-1; dx<=+1; ++dx) {
        for (plint dy=-1; dy<=+1; ++dy) {
            for (plint dz=-1; dz<=+1; ++dz) {
                if (!(dx==0 && dy==0 && dz==0)) {
                    Dot3D neigh = pos+offset+Dot3D(dx,dy,dz);
                    int typeOfNeighbor = voxels.get(pos.x+dx,pos.y+dy,pos.z+dz);
                    if (typeOfNeighbor!=voxelFlag::undetermined) {
                        T distance;
                        plint whichTriangle;
                        Array<T,3> p1(pos_.x,pos_.y,pos_.z);
                        Array<T,3> p2(neigh.x,neigh.y,neigh.z);
                        global::counter("voxelizer-debug").increment(1);
                        bool crossed = checkIfFacetsCrossed (hashContainer, p1, p2, distance, whichTriangle);
                        global::counter("voxelizer-debug").reset();
                        std::cout << "Neighbor ("
                                  << dx << "," << dy << "," << dz
                                  << "); is "
                                  << (voxelFlag::insideFlag(typeOfNeighbor) ? "inside" : "outside");
                        if (crossed) {
                            triangles.insert(whichTriangle);
                            std::cout 
                                  << " inters. at distance " << distance
                                  << " with triangle " << whichTriangle << std::endl;
                        }
                        else {
                            std::cout << " no inters." << std::endl;
                        }
                    }
                }
            }
        }
    }
    std::set<plint>::iterator it = triangles.begin();
    for (; it!=triangles.end(); ++it) {
        std::cout << "Triangle " << *it << " [" << std::flush;
        Array<T,3> p0 = mesh.getVertex(*it, 0);
        Array<T,3> p1 = mesh.getVertex(*it, 1);
        Array<T,3> p2 = mesh.getVertex(*it, 2);
        std::cout << p0[0] << " " << p1[0] << " " << p2[0] << " " << p0[0] << "], ["
                  << p0[1] << " " << p1[1] << " " << p2[1] << " " << p0[1] << "], ["
                  << p0[2] << " " << p1[2] << " " << p2[2] << " " << p0[2] << "]" << std::endl;
    }
}

template<typename T>
bool VoxelizeMeshFunctional3D<T>::voxelizeFromNeighbor(ScalarField3D<int> const& voxels, AtomicContainerBlock3D& hashContainer,
	const Dot3D& pos, const Dot3D& neighbor, int& voxelType, const int& verificationLevel)
{
    Dot3D offset = voxels.getLocation();
    int typeOfNeighbor = voxels.get(neighbor.x,neighbor.y,neighbor.z);
    if (typeOfNeighbor==voxelFlag::undetermined) {
        return true;
    }
    // If there is no verification and the voxel has already been voxelized,
    //   it is not being re-voxelized here.
    if (verificationLevel==0) {
        if (voxelType!=voxelFlag::undetermined) {
            return true;
        }
    }
    Dot3D pos_ = pos+offset;
    Dot3D neighbor_ = neighbor+offset;
    Array<T,3> point1((T)pos_.x, (T)pos_.y, (T)pos_.z);
    Array<T,3> point2((T)neighbor_.x, (T)neighbor_.y, (T)neighbor_.z);
    int newVoxelType = voxelFlag::undetermined;
    T distance1, distance2, distance3, distance4;
    bool isBehind1, isBehind2;
    plint whichTriangle1, whichTriangle2;
    if (checkIfFacetsCrossed(hashContainer, point1, point2, distance1, whichTriangle1)) {
        newVoxelType = voxelFlag::invert(typeOfNeighbor);
        // Additional consistency checks only at the ultimate level of verification.
        if (verificationLevel==2) {
            PLB_ASSERT( distance1 < std::sqrt((T)3)+(T)0.0001 );
#ifdef PLB_DEBUG
            bool ok = checkIfFacetsCrossed(hashContainer, point2, point1, distance2, whichTriangle2);
#else
            (void) checkIfFacetsCrossed(hashContainer, point2, point1, distance2, whichTriangle2);
#endif
            PLB_ASSERT( ok );
            PLB_ASSERT( distance2 < std::sqrt((T)3)+(T)0.0001 );

#ifdef PLB_DEBUG
            bool ok1 = distanceToSurface( hashContainer, point1, distance3, isBehind1);
#else
            (void) distanceToSurface( hashContainer, point1, distance3, isBehind1);
#endif

            PLB_ASSERT( ok1 );
            PLB_ASSERT( distance1 < std::sqrt((T)3)+(T)0.0001 );
            // Attention: At this moment, the following consistency check fails sometimes,
            //   god knows why. It might be that there is a bug in the method
            //   mesh.distanceToSurface.
            PLB_ASSERT( (voxelFlag::insideFlag(newVoxelType) && isBehind1) ||
                        (voxelFlag::outsideFlag(newVoxelType) && !isBehind1) );

#ifdef PLB_DEBUG
            bool ok2 = distanceToSurface( hashContainer, point2, distance4, isBehind2);
#else
            (void) distanceToSurface( hashContainer, point2, distance4, isBehind2);
#endif
            PLB_ASSERT( ok2 );
            PLB_ASSERT( distance2 < std::sqrt((T)3)+(T)0.0001 );
            PLB_ASSERT ( (voxelFlag::insideFlag(typeOfNeighbor) && isBehind2) ||
                         (voxelFlag::outsideFlag(typeOfNeighbor) && !isBehind2) );
        }
    }
    else {
        newVoxelType = typeOfNeighbor;
    }
    int oldVoxelType = voxelType;
    voxelType = newVoxelType;
    if (oldVoxelType == voxelFlag::undetermined) {
        return true;
    }
    else {
        return oldVoxelType == newVoxelType;
    }
}

template<typename T>
void VoxelizeMeshFunctional3D<T>::processGenericBlocks (Box3D domain, const std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION( blocks.size()==2 );
    ScalarField3D<int>* voxels =	dynamic_cast<ScalarField3D<int>*>(blocks[0]);
    PLB_ASSERT( voxels );
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[1]);
    PLB_ASSERT( container );

	#ifdef PLB_DEBUG
		bool main = false;
		main = global::mpi().isMainProcessor();
		if(main){std::cout << "[DEBUG] VoxelizeMeshFunctional3D<T>::processGenericBlocks" << std::endl;}
	#endif

    // Return if this block is already voxelized.
    if (voxels->getFlag()) {
        return;
    }

    Array<plint,2> xRange, yRange, zRange;
    if (!createVoxelizationRange(domain, *voxels, xRange, yRange, zRange)) {
        // If no seed has been found in the envelope, just return and wait
        //   for the next round.
        return;
    }

	plint minX = xRange[1]>xRange[0] ? xRange[0] : xRange[1];
	minX++;
	plint maxX = xRange[1]>xRange[0] ? xRange[1] : xRange[0];
	maxX--;
	plint minY = yRange[1]>yRange[0] ? yRange[0] : yRange[1];
	minY++;
	plint maxY = yRange[1]>yRange[0] ? yRange[1] : yRange[0];
	maxY--;
	plint minZ = zRange[1]>zRange[0] ? zRange[0] : zRange[1];
	minZ++;
	plint maxZ = zRange[1]>zRange[0] ? zRange[1] : zRange[0];
	maxZ--;

	#ifdef PLB_DEBUG
		if(main){std::cout << "[DEBUG] Finding rotten Voxels" << std::endl;}
	#endif

	#ifdef PLB_MPI_PARALLEL
		const int nproc = global::mpi().getSize();
		const int rank = global::mpi().getRank();
		std::vector<Box3D> mpiDomains = global::mpiData().splitDomains(domain); //Overwrite the given domain
		domain = mpiDomains[rank];
		minX = domain.x0;
		maxX = domain.x1;
		minY = domain.y0;
		maxY = domain.y1;
		minZ = domain.z0;
		maxZ = domain.z1;
	#endif

	std::vector<Dot3D> undeterminedVoxels;

	const plint nBlocks = (maxX-minX)*(maxY-minY)*(maxZ-minZ);
	undeterminedVoxels.resize(nBlocks);

	for(plint iX = minX; iX < maxX; iX++){
		for (plint iY = minY; iY < maxY; iY++) {
			for (plint iZ = minZ; iZ < maxZ; iZ++) {
				if(voxels->get(iX, iY, iZ) == voxelFlag::undetermined){
					Dot3D dot(iX, iY, iZ);
					undeterminedVoxels.push_back(dot);
				}
			}
		}
	}

	const plint nVoxels = undeterminedVoxels.size();
	std::vector<std::vector<Dot3D> > voxelRepair(nVoxels);

	#ifdef PLB_DEBUG
		if(main){std::cout << "[DEBUG] Finding healthy neighbours for "<<nVoxels<<" voxels." << std::endl;}
	#endif

	int nx = voxels->getNx(); int ny = voxels->getNy(); int nz = voxels->getNz();
	for(int n=0; n<=nVoxels; n++){
		Dot3D pos = undeterminedVoxels[n];
		int x = pos.x; int y = pos.y; int z = pos.z;
		if((x >= nx) || (x <= 0) || (y >= ny) || (y <= 0) || (z >= nz) || (z <= 0)){
			continue;
		}
		int voxelType = voxels->get(pos.x,pos.y,pos.z);
		if(voxelType == voxelFlag::undetermined){
			std::vector<Dot3D> neighbours;
			for (plint dx=-1; dx<=+1; ++dx) {
				for (plint dy=-1; dy<=+1; ++dy) {
					for (plint dz=-1; dz<=+1; ++dz) {
						if(dx==0 && dy==0 && dz==0){ continue;}
						else{
							x = pos.x+dx; y = pos.y+dy; z = pos.z+dz;
							if((x >= nx) || (x <= 0) || (y >= ny) || (y <= 0) || (z >= nz) || (z <= 0)){
								if(main){std::cout << "Neighbour Voxels.get(x,y,z) Offender = "<< x << ", " << y << ", " << z << std::endl;}
								continue;
							}
							else{ if(voxels->get(x, y, z)!=voxelFlag::undetermined){neighbours.push_back(Dot3D(x, y, z));} }
						}
					}
				}
			}
			voxelRepair[n] = neighbours;
		}
	}

	#ifdef PLB_DEBUG
		if(main){std::cout << "[DEBUG] Fixing "<<nVoxels<<" voxels." << std::endl;}
	#endif

	int verificationLevel = 0;

	for(int n=0; n<=nVoxels; n++){
		if(n>voxelRepair.size()-1){ break; }
		bool offender = false;
		Dot3D pos = undeterminedVoxels[n];
		int x = pos.x; int y = pos.y; int z = pos.z;
		if((x >= nx) || (x <= 0) || (y >= ny) || (y <= 0) || (z >= nz) || (z <= 0)){
			if(main){std::cout << "Fix Voxels.get(x,y,z) Offender = "<< x << ", " << y << ", " << z << std::endl;}
			offender = true;
		}
		if(!offender){
			int voxelType = voxels->get(pos.x,pos.y,pos.z);
			std::vector<Dot3D> neighbours = voxelRepair[n];
			plint nb = neighbours.size();
			for(int i = 0; i<nb; i++)
			{
				bool ok = voxelizeFromNeighbor (*voxels, *container, pos, neighbours[i], voxelType, verificationLevel);
				if (!ok) { printOffender(*voxels, *container, pos);}
				else{ break; }
			}
			voxels->get(pos.x, pos.y, pos.z) = voxelType;
		}
	}

	#ifdef PLB_MPI_PARALLEL
		#ifdef PLB_DEBUG
			if(main){std::cout << "[DEBUG] Sharing voxel data across mpi processes" << std::endl;}
		#endif
		// Merge results in the voxels ScalarField
		for(int id = 0; id < nproc; id++){
			if(id == rank){ global::mpiData().sendScalarField3D(*voxels,domain); }
			else{global::mpiData().receiveScalarField3D(*voxels,mpiDomains[id],rank);}
		}
	#endif

    // Indicate that this atomic-block has been voxelized.
	voxels->setFlag(true);
	#ifdef PLB_DEBUG
		if(main){std::cout << "[DEBUG] DONE VoxelizeMeshFunctional3D<T>::processGenericBlocks" << std::endl;}
	#endif
}

template<typename T>
VoxelizeMeshFunctional3D<T>* VoxelizeMeshFunctional3D<T>::clone() const {
    return new VoxelizeMeshFunctional3D<T>(*this);
}

template<typename T>
void VoxelizeMeshFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Voxels
    modified[1] = modif::nothing; // Hash Container
}

template<typename T>
BlockDomain::DomainT VoxelizeMeshFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** DetectBorderLineFunctional3D ************************************* */

template<typename T>
void detectBorderLine( MultiScalarField3D<T>& voxelMatrix, Box3D const& domain, plint borderWidth ){
    applyProcessingFunctional( new DetectBorderLineFunctional3D<T>(borderWidth), domain, voxelMatrix );
}

template<typename T>
DetectBorderLineFunctional3D<T>::DetectBorderLineFunctional3D(plint borderWidth_): borderWidth(borderWidth_){ }

template<typename T>
void DetectBorderLineFunctional3D<T>::process (Box3D domain, ScalarField3D<T>& voxels )
{
	#ifdef PLB_DEBUG
		std::cout<< "[DEBUG] DetectBorderLineFunctional3D<T>::process"<< std::endl;
	#endif
	Box3D voxelBox = voxels.getBoundingBox();
	plint size = voxels.getSize();

	const plint minX = domain.x1>domain.x0 ?  domain.x0 - borderWidth :  domain.x1 - borderWidth;
	const plint maxX = domain.x1>domain.x0 ?  domain.x1 + borderWidth :  domain.x0 + borderWidth;

	const plint minY = domain.y1>domain.y0 ?  domain.y0 - borderWidth :  domain.y1 - borderWidth;
	const plint maxY = domain.y1>domain.y0 ?  domain.y1 + borderWidth :  domain.y0 + borderWidth;

	const plint minZ = domain.z1>domain.z0 ?  domain.z0 - borderWidth :  domain.z1 - borderWidth;
	const plint maxZ = domain.z1>domain.z0 ?  domain.z1 + borderWidth :  domain.z0 + borderWidth;

	std::vector<Dot3D> borderVoxels;

	for(plint iX=minX; iX<maxX; iX++){
		for (plint iY = minY; iY <= maxY; iY++) {
			for (plint iZ = minZ; iZ <= maxZ; iZ++) {
				int voxelType = voxels.get(iX, iY, iZ);
				if(voxelFlag::outsideFlag(voxelType) || voxelFlag::insideFlag(voxelType)){
					Dot3D dot(iX, iY, iZ);
					borderVoxels.push_back(dot);
				}
			}
		}
	}

	size = borderVoxels.size();

	for(int n=0; n<=size; n++){
		Dot3D dot = borderVoxels[n];
		for (plint dx=-borderWidth; dx<=borderWidth; ++dx){
			for (plint dy=-borderWidth; dy<=borderWidth; ++dy){
				for (plint dz=-borderWidth; dz<=borderWidth; ++dz){
					if(!(dx==0 && dy==0 && dz==0)) {
						plint nextX = dot.x + dx;
						plint nextY = dot.y + dy;
						plint nextZ = dot.z + dz;
						Dot3D neighbour(nextX,nextY,nextZ);
						if(contained(neighbour,voxelBox)){
							int dotFlag = voxels.get(dot.x, dot.y, dot.z);
							int neighbourFlag = voxels.get(neighbour.x, neighbour.y, neighbour.z);
							if(voxelFlag::outsideFlag(dotFlag) && voxelFlag::insideFlag(neighbourFlag)){
								voxels.get(dot.x,dot.y,dot.z) = voxelFlag::outerBorder;}
							if(voxelFlag::insideFlag(dotFlag) && voxelFlag::outsideFlag(neighbourFlag)){
								voxels.get(dot.x,dot.y,dot.z) = voxelFlag::innerBorder;}
						}
					}
				}
			}
		}
	}

	#ifdef PLB_DEBUG
		std::cout<< "[DEBUG] DONE DetectBorderLineFunctional3D<T>::process"<< std::endl;
	#endif
}

template<typename T>
DetectBorderLineFunctional3D<T>* DetectBorderLineFunctional3D<T>::clone() const {
    return new DetectBorderLineFunctional3D<T>(*this);
}

template<typename T>
void DetectBorderLineFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT DetectBorderLineFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

} // namespace plb

#endif  // VOXELIZER_HH
