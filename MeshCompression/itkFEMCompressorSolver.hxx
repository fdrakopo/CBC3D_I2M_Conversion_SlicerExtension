/*
 * itkFEMCompressorSolver.hxx
 *
 */

#ifndef ITKFEMCOMPRESSORSOLVER_HXX_
#define ITKFEMCOMPRESSORSOLVER_HXX_


#include "itkFEMCompressorSolver.h"
#include "itkFEMLoadNoisyLandmark.h"

// vtk
#include "vtkUnstructuredGrid.h"
#include "vtkIdList.h"
#include "vtkUnstructuredGridWriter.h"

namespace itk
{
namespace fem
{

template <class TPixel,unsigned int VDimension>
CompressorSolver<TPixel,VDimension>
::CompressorSolver()
 {
	this->m_ForceIndex = 0;
	this->m_SolutionIndex = 0;
	this->m_MeshStiffnessMatrixIndex = 1;
	this->m_LandmarkStiffnessMatrixIndex = 2;
	this->m_StiffnessMatrixIndex = 0;

	this->m_SurfaceVertexToLabelMapContainer = NULL;
	this->m_SurfaceVertexToNormalMapContainer = NULL;
	this->m_SurfaceVertexToCorrespondenceMapContainer = NULL;
	this->m_SurfaceVertexToNeighborVerticesMapContainer = NULL;
	this->m_SurfaceTrianglesContainer = NULL;
	this->m_TargetPointSet = NULL;
	this->m_SurfaceVertexToNeighborTargetPointsMapContainer = NULL;
	this->m_TetrahedraDataContainer = NULL;

	this->m_NumIterations = 5;
	this->m_Lambda = 1.0; // (0,1]
	this->m_ScaleFactorForEdges = 0.25; // (0,1]
	this->m_ScaleFactorForVolumes = 0.25; // (0,1]
	this->m_ScaleFactorForSearchRadius = 2.0; // > 1

 }

template <class TPixel,unsigned int VDimension>
CompressorSolver<TPixel,VDimension>
::~CompressorSolver()
{
}

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::GenerateData()
 {
	// (0,1]
	if(this->GetScaleFactorForVolumes() <= 0 || this->GetScaleFactorForVolumes() > 1)
		itkExceptionMacro("Invalid value for Scale Factor For Volumes : " << this->GetScaleFactorForVolumes() << " .Should be between (0,1].");

	// (0,1]
	if(this->GetScaleFactorForEdges() <= 0 || this->GetScaleFactorForEdges() > 1)
		itkExceptionMacro("Invalid value for Scale Factor For Edges : " << this->GetScaleFactorForEdges() << " .Should be between (0,1].");

	// > 1
	if(this->GetScaleFactorForSearchRadius() < 1)
		itkExceptionMacro("Invalid value for Scale Factor For Search Radius : " << this->GetScaleFactorForSearchRadius() << " .Should be > 1.");

	// (0,1]
	if(this->GetLambda() <= 0 || this->GetLambda() > 1.0)
		itkExceptionMacro("Invalid value for Lamdba : " << this->GetLambda() << " .Should be between (0,1].");

	//std::cout << "Scale Factor For Edges : " << this->m_ScaleFactorForEdges << std::endl;
	//std::cout << "Scale Factor For Volumes : " << this->m_ScaleFactorForVolumes << std::endl;
	//std::cout << "Scale Factor For Search Radious : " << this->m_ScaleFactorForSearchRadius << std::endl;

	// Initialize matrix, vector, solution
	this->Initialization();

	// Compute for each source point the neighbor target points
	this->ComputeNeighborTargetPoints(this->GetScaleFactorForSearchRadius());

	// Add each surface vertex to it's own neighbor set.
	// We avoid insertion and removal in the compressor loop.
	this->AddSurfaceVertexToItsOwnNeighborsSet();

	std::stringstream iter;
	std::string fname;
	std::string meshName = "mesh_iter_";

	unsigned int numIter = this->GetNumIterations();

	// The Compressor loop
	for(unsigned int i = 0; i < numIter; i++)
	{
		std::cout << "\nIteration : " << i+1 << " / " << numIter << std::endl;

		this->AssembleMeshStiffnessMatrix();
		this->ComputeLandmarksTensor();
		this->AssembleLandmarkStiffnessMatrix();
		this->AssembleGlobalMatrixFromLandmarksAndMeshMatrices();
		this->ComputeNormals();
		this->NormalizeNormals();
		this->ComputeCorrespondence();
		this->UpdateLandmarks();
		this->AssembleF();
		// this->PrintForces("forces.txt");
		//std::cout << "Solving System..." << std::endl;
		this->GetLinearSystemWrapper()->Solve();
		// this->PrintSolution("solution.txt");
		this->UpdateFEMObjectCoordinates();

		// Check Tetrahedra Volumes
		if(!CheckTetrahedraVolumes(this->m_FEMObject->GetElementContainer()))
		{
			// If found inverted tetrahedra we restore the previous mesh ccordinates and exit the loop
			RestorePreviousFEMObjectCoordinates();
			// if(!CheckTetrahedraVolumes(this->m_FEMObject->GetElementContainer()))
			// 	itkExceptionMacro("Restore Previous FEMObject Coordinates Failed!");
			break;
		}

//		// Save Tets
//		iter.str("");
//		iter << i+1;
//		fname = meshName + iter.str() + ".vtk";
//		SaveTetrahedraVTK(this->m_FEMObject->GetNodeContainer(),this->m_FEMObject->GetElementContainer(),this->GetTetrahedraDataContainer(),fname);
	}

	// Remove each surface vertex from it's own neighbor set.
	this->RemoveSurfaceVertexFromItsOwnNeighborsSet();

	// The copy from input to output takes to much time especially for large meshes
	// FEMObjectType* pOutputFEMObject = this->GetOutput();
	// pOutputFEMObject->DeepCopy(this->GetInput());
	// pOutputFEMObject->FinalizeMesh();
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel, VDimension>
::ComputeNeighborTargetPoints(const double fac)
 {
	std::cout << "Computing Neighbor Target Points..." <<  std::endl;

	if(!this->m_SurfaceVertexToLabelMapContainer)
		itkExceptionMacro("Surface Vertex To Label Map Container is NULL!");
	if(!this->m_TargetPointSet || !this->m_TargetPointSet->GetPoints())
		itkExceptionMacro("Target Points is NULL!");
	unsigned int numTargetPoints = this->m_TargetPointSet->GetPoints()->Size();
	if(numTargetPoints < 1)
		itkExceptionMacro("Invalid number of target points.");

	this->m_SurfaceVertexToNeighborTargetPointsMapContainer = SurfaceVertexToNeighborTargetPointsMapContainerType::New();
	this->m_SurfaceVertexToNeighborTargetPointsMapContainer->Initialize();

	typename SurfaceVertexToLabelMapContainerType::Iterator sourceIt;
	unsigned int id = 0;
	VectorType sPos,sNeighborPos;
	VectorType sDist;
	typename TargetPointSetType::PointType tPt;
	PixelTypeSet sLabels,tLabels;

	// Loop over the surface vertices
	for(sourceIt = this->m_SurfaceVertexToLabelMapContainer->Begin(); sourceIt != this->m_SurfaceVertexToLabelMapContainer->End();  sourceIt++)
	{
		id = sourceIt->Index();
		sPos = this->m_FEMObject->GetNode(id)->GetCoordinates();
		sLabels = this->m_SurfaceVertexToLabelMapContainer->ElementAt(id);

		// Get the source neighbors of the vertex.
		std::set<unsigned int>* pSourceNeighbors = this->m_SurfaceVertexToNeighborVerticesMapContainer->ElementAt(id);
		if(!pSourceNeighbors)
			itkExceptionMacro("The neighbors vertices vector is NULL!");
		if(pSourceNeighbors->size() < 2)
			itkExceptionMacro("The size of the neighbor set is less than 2.");

		std::set<unsigned int>::iterator neighborIt;
		double searchRadius = 0;
		unsigned int sId = 0;

		// Loop through all the source neighbor vertices
		for(neighborIt = pSourceNeighbors->begin(); neighborIt != pSourceNeighbors->end(); neighborIt++)
		{
			sId = *neighborIt;
			sNeighborPos = this->m_FEMObject->GetNode(sId)->GetCoordinates();
			sDist = sNeighborPos - sPos;
			searchRadius += sDist.magnitude();
		}

		// The search radious is magnified by the factor
		searchRadius /= pSourceNeighbors->size();
		searchRadius *= fac;

		// Construct the box with half sides equal to radius
		VectorType minBoxCoord = sPos;
		VectorType maxBoxCoord = sPos;
		for(unsigned int i=0; i<3; i++)
		{
			minBoxCoord[i] -= searchRadius;
			maxBoxCoord[i] += searchRadius;
		}

		std::vector<unsigned int>* pTargetNeighbors = new std::vector<unsigned int>();

		bool bInsideBox = false;

		// Loop over all target points
		for(unsigned int i=0; i < numTargetPoints; i++)
		{
			tPt = this->m_TargetPointSet->GetPoint(i);
			tLabels = this->m_TargetPointSet->GetPointData()->ElementAt(i);

			// Check if inside the search bounding box
			if(!IsInsideBox(tPt,minBoxCoord,maxBoxCoord))
				continue;
			else
				bInsideBox = true;

			// Check only points with same label sets
			if(!SameLabelSets(sLabels,tLabels))
				continue;

			pTargetNeighbors->push_back(i);
		}

		// Search again with a partial match of the image labels (only for the case of 2 labels)
		if(bInsideBox && pTargetNeighbors->size() == 0 && sLabels.size() == 2)
		{
			for(unsigned int i=0; i < numTargetPoints; i++)
			{
				tPt = this->m_TargetPointSet->GetPoint(i);
				tLabels = this->m_TargetPointSet->GetPointData()->ElementAt(i);

				if(tLabels.size() != 2)
					continue;

				// Check if inside the search bounding box
				if(!IsInsideBox(tPt,minBoxCoord,maxBoxCoord))
					continue;

				// Check only if partially the labels are same
				if(!PartiallySameLabelSets(sLabels,tLabels))
					continue;

				pTargetNeighbors->push_back(i);
			}
		}

		if(pTargetNeighbors->size() == 0)
		{
			std::cout << "Warning : Target Neighbors size is zero, check search radious or labels!" << std::endl;
			std::cout << "source Labels : " << std::endl;
			typename PixelTypeSet::iterator it;
			for(it = sLabels.begin(); it != sLabels.end(); ++it)
				std::cout << ' ' << *it;
			std::cout << std::endl;
		}

//		// Print
//		if(sLabels.first != 0 && sLabels.second != 0)
//		{
//			PrintVertex(sPos,"sourcePoint.vtk");
//			PrintNeighborTargetVertices(pTargetNeighbors,this->m_TargetPointSet->GetPoints(),"targetNeighbors.vtk");
//		}

		// insert the pair to the map
		this->m_SurfaceVertexToNeighborTargetPointsMapContainer->InsertElement(id,pTargetNeighbors);
	}
	//std::cout << "Size : " << this->m_SurfaceVertexToNeighborTargetPointsMapContainer->Size() << std::endl;
 }


template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel, VDimension>
::PrintNeighborTargetVertices(std::vector<unsigned int>* pVector,typename TargetPointSetType::PointsContainer* pPointsContainer,std::string fName)
{
	std::cout << "Saving Vertices at  : " << fName << std::endl;

	if(!pVector || !pPointsContainer || fName.empty())
	{
		std::cout << "Invalid Input!" << std::endl;
		return;
	}

	vtkUnstructuredGrid * pUG = vtkUnstructuredGrid::New();
	vtkPoints* pPoints = vtkPoints::New();

	// Set coordinates
	int numPoints = pPointsContainer->Size();

	typename TargetPointSetType::PointType pt;

    for( int i=0; i < numPoints; i++ )
    {
    	pt = pPointsContainer->ElementAt(i);
    	pPoints->InsertPoint(i,pt[0],pt[1],pt[2]);
    }

    // Set the points to Unstructured Grid
	pUG->SetPoints(pPoints);
	pPoints->Delete();

	vtkIdList* ptIds = vtkIdList::New();
	long int id;

    int numVertices = pVector->size();

	for(int i = 0; i < numVertices; i++)
	{
		id = pVector->at(i);
		ptIds->InsertId(0,id);

		// Insert the cell to Unstructured Grid
		pUG->InsertNextCell(1,ptIds);
	}
	ptIds->Delete();

	vtkUnstructuredGridWriter* pWriter = vtkUnstructuredGridWriter::New();
	pWriter->SetFileTypeToASCII();
	pWriter->SetFileName(fName.c_str());
	pWriter->SetInputData(pUG);
	pWriter->Write();
	pWriter->Update();
	pWriter->Delete();
}


template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel, VDimension>
::UpdateFEMObjectCoordinates()
 {
	//std::cout << "Updating FEMObject Node Coordinates..." << std::endl;

	int numNodes = this->m_FEMObject->GetNumberOfNodes();
	typedef Element::Node NodeType;
	VectorType pt(FEMDimension);

	for( int i = 0; i < numNodes; i++ )
	{
		NodeType::Pointer node = this->m_FEMObject->GetNode(i);
		for( unsigned int j = 0; j < FEMDimension; j++ )
			pt[j] = node->GetCoordinates()[j] + this->GetLinearSystemWrapper()->GetSolutionValue(node->GetDegreeOfFreedom(j));

		node->SetCoordinates(pt);
	}
}

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel, VDimension>
::RestorePreviousFEMObjectCoordinates()
 {
	//std::cout << "Restoring Previous FEMObject Coordinates..." << std::endl;

	int numNodes = this->m_FEMObject->GetNumberOfNodes();
	typedef Element::Node NodeType;
	VectorType pt(FEMDimension);

	for( int i = 0; i < numNodes; i++ )
	{
		NodeType::Pointer node = this->m_FEMObject->GetNode(i);
		for( unsigned int j = 0; j < FEMDimension; j++ )
			pt[j] = node->GetCoordinates()[j] - this->GetLinearSystemWrapper()->GetSolutionValue(node->GetDegreeOfFreedom(j));

		node->SetCoordinates(pt);
	}
}


template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel, VDimension>
::NormalizeNormals()
 {
	//std::cout << "Normalizing Normals..." << std::endl;
	if(!this->m_SurfaceVertexToNormalMapContainer)
		itkExceptionMacro("Surface Triangles Container is NULL!");

	unsigned int vId;
	VectorType normal;

	typename SurfaceVertexToNormalMapContainerType::Iterator iter;
	for(iter = this->m_SurfaceVertexToNormalMapContainer->Begin(); iter != this->m_SurfaceVertexToNormalMapContainer->End(); iter++)
	{
		vId = iter->Index();
		normal = iter->Value();
		normal.normalize();
		this->m_SurfaceVertexToNormalMapContainer->ElementAt(vId) = normal;
	}
 }


template <class TPixel,unsigned int VDimension>
typename itk::fem::CompressorSolver<TPixel,VDimension>::VectorType
CompressorSolver<TPixel,VDimension>
::ComputeTriangleNormal(VectorTriIdType vTriId)
 {
	VectorType normal;
	normal.fill(0);

	// Get the coordinates of the vertices
	VectorType pos0 = this->m_FEMObject->GetNode(vTriId[0])->GetCoordinates();
	VectorType pos1 = this->m_FEMObject->GetNode(vTriId[1])->GetCoordinates();
	VectorType pos2 = this->m_FEMObject->GetNode(vTriId[2])->GetCoordinates();

	VectorType vec01 = pos1 - pos0;
	VectorType vec02 = pos2 - pos0;

	normal = vnl_cross_3d(vec01,vec02);
	//std::cout << normal << std::endl;
	normal.normalize();
	//std::cout << normal << std::endl;
	return normal;
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::ComputeNormals()
 {
	//std::cout << "Computing Normals..." << std::endl;
	if(!this->m_SurfaceVertexToNormalMapContainer)
		this->m_SurfaceVertexToNormalMapContainer = SurfaceVertexToNormalMapContainerType::New();

	// Tell the container to release any memory it may have allocated and return itself to its initial state.
	this->m_SurfaceVertexToNormalMapContainer->Initialize();

	if(!this->m_SurfaceTrianglesContainer)
		itkExceptionMacro("Surface Triangles Container is NULL!");

	typename SurfaceTrianglesContainerType::Iterator iter;
	VectorTriIdType vTriId;
	VectorType normal;
	std::pair<unsigned int,VectorType> vertexToNormal_pair;
	typename SurfaceVertexToNormalMapContainerType::iterator findIt;

	// Loop over all surface triangles
	for(iter = this->m_SurfaceTrianglesContainer->Begin(); iter != this->m_SurfaceTrianglesContainer->End(); iter++)
	{
		// The vertices ids of the surface triangle
		vTriId = iter.Value();

		// Compute the normal
		normal = this->ComputeTriangleNormal(vTriId);
		vertexToNormal_pair.second = normal;

		for(unsigned int i = 0; i < FEMDimension; i++)
		{
			vertexToNormal_pair.first = vTriId[i];
			findIt = this->m_SurfaceVertexToNormalMapContainer->CastToSTLContainer().find(vTriId[i]);

			// The vertex id exist so we add the normal to the existing value of the map
			if(findIt != this->m_SurfaceVertexToNormalMapContainer->CastToSTLContainer().end())
			{
				//std::cout << findIt->second << std::endl;
				this->m_SurfaceVertexToNormalMapContainer->CastToSTLContainer()[vTriId[i]] = findIt->second + vertexToNormal_pair.second;
				//std::cout << this->m_SurfaceVertexToNormalMapContainer->CastToSTLContainer()[vTriId[i]] << std::endl;
			}
			// The vertex id doesn't exist
			else
			{
				this->m_SurfaceVertexToNormalMapContainer->CastToSTLContainer().insert(vertexToNormal_pair);
			}
		}
	}
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::UpdateLandmarks()
 {
	//std::cout << "Updating Landmarks..." << std::endl;

	if(!this->m_FEMObject->GetLoadContainer())
		itkExceptionMacro("Missing load container!");

	// Check
	if(this->m_FEMObject->GetLoadContainer()->Size() != this->m_SurfaceVertexToCorrespondenceMapContainer->Size())
		itkExceptionMacro("Incompatible container sizes between Load Container and Surface Vertex To Correspondence Map Container!");

	typename SurfaceVertexToCorrespondenceMapContainerType::Iterator iter;
	VectorType correspondence;
	unsigned int count=0;
	unsigned int id=0;
	VectorType pos;
	Load::Pointer pLoad;
	LoadType* pLoadNoisyLandmark;

	for(iter = this->m_SurfaceVertexToCorrespondenceMapContainer->Begin(); iter != this->m_SurfaceVertexToCorrespondenceMapContainer->End(); iter++)
	{
		id = iter->Index();
		correspondence = iter->Value();

		// The current pos of the mesh surface node
		pos = this->m_FEMObject->GetNode(id)->GetCoordinates();

		pLoad = this->m_FEMObject->GetLoadContainer()->ElementAt(count);
		pLoadNoisyLandmark = dynamic_cast<LoadType*>( pLoad.GetPointer() );
		pLoadNoisyLandmark->SetSource(pos);
		pLoadNoisyLandmark->SetRealDisplacement(correspondence);
		count++;
	}
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::ComputeCorrespondence()
 {
	//std::cout << "Computing Correspondences..." << std::endl;

	if(!this->m_SurfaceVertexToCorrespondenceMapContainer)
		this->m_SurfaceVertexToCorrespondenceMapContainer = SurfaceVertexToCorrespondenceMapContainerType::New();

	// Tell the container to release any memory it may have allocated and return itself to its initial state.
	this->m_SurfaceVertexToCorrespondenceMapContainer->Initialize();

	if(!this->m_SurfaceVertexToLabelMapContainer)
		itkExceptionMacro("Surface Vertex To Label Map Container is NULL!");
	if(!this->m_SurfaceVertexToNormalMapContainer)
		itkExceptionMacro("Surface Vertex To Normal Map Container is NULL!");
	if(!this->m_SurfaceVertexToNeighborVerticesMapContainer)
		itkExceptionMacro("Surface Vertex To Neighbor Vertices Map Container is NULL!");
	if(!this->m_SurfaceVertexToNeighborTargetPointsMapContainer)
		itkExceptionMacro("Surface Vertex To Neighbor Target Points is NULL!");
	if(!this->m_TargetPointSet)
		itkExceptionMacro("Target Point Set is NULL!");

	typename SurfaceVertexToLabelMapContainerType::Iterator vertexToLabelIt;

	VectorType sPos;
	VectorType tPos(FEMDimension);
	typename TargetPointSetType::PointType tPt;
	VectorType closest_tPos(FEMDimension);
	VectorType sum_tPos(FEMDimension),sum_tPos_proj(FEMDimension);

	VectorType distPos,distPos_proj;
	double dist=0;
	VectorType normal;
	VectorType zeroPos(FEMDimension);
	zeroPos.fill(0);

	unsigned int id = 0;
	unsigned int sId = 0;

	int numSourceNeighbors = 0;
	const double min_dist_init = std::numeric_limits<double>::max();

	// Loop over the surface vertices
	for(vertexToLabelIt = this->m_SurfaceVertexToLabelMapContainer->Begin(); vertexToLabelIt != this->m_SurfaceVertexToLabelMapContainer->End();  vertexToLabelIt++)
	{
		id = vertexToLabelIt->Index();

		// Get the source neighbors of the vertex.
		// We have already INSERTED the surface vertex to it's own Neighbor set.
		std::set<unsigned int>* pSourceNeighbors = this->m_SurfaceVertexToNeighborVerticesMapContainer->ElementAt(id);
		if(!pSourceNeighbors)
			itkExceptionMacro("The Neighbor Source Points vector is NULL!");
		if(pSourceNeighbors->size() < 3)
			itkExceptionMacro("The size of the neighbor set is less than 3!");

		// This is including the current surface vertex
		numSourceNeighbors = pSourceNeighbors->size();

		sum_tPos = zeroPos;
		std::set<unsigned int>::iterator neighborSourceIt;

		// For each source neighbor vertex found the closest target with the same label
		for(neighborSourceIt = pSourceNeighbors->begin(); neighborSourceIt != pSourceNeighbors->end(); neighborSourceIt++)
		{
			sId = *neighborSourceIt;
			// if(sId==3182 && id==2833)
			// int d = 84;

			sPos = this->m_FEMObject->GetNode(sId)->GetCoordinates();
			//sLabels = this->m_SurfaceVertexToLabelMapContainer->ElementAt(sId);

			double min_dist = min_dist_init;
			bool bFound = false;

			std::vector<unsigned int>::iterator neighborTargetIt;
			std::vector<unsigned int>* pTargetNeighbors = this->m_SurfaceVertexToNeighborTargetPointsMapContainer->ElementAt(sId);
			if(!pTargetNeighbors)
				itkExceptionMacro("Neighbor Target Points vector is NULL!");

			// Loop over all the neighbor target points of the source point
			for(neighborTargetIt = pTargetNeighbors->begin(); neighborTargetIt != pTargetNeighbors->end(); neighborTargetIt++)
			{
				tPt = this->m_TargetPointSet->GetPoint(*neighborTargetIt);
				//tLabels = this->m_TargetPointSet->GetPointData()->ElementAt(*neighborTargetIt);

				//				// Here the source and the target point have the same labels
				//				if(!SameLabels(sLabels,tLabels))
				//					itkExceptionMacro("Source and target point have different labels!");

				for(unsigned int j=0; j < FEMDimension; j++)
					tPos[j] = tPt[j];

				distPos = tPos - sPos;
				//std::cout << distPos << std::endl;
				dist = distPos.two_norm();
				if(dist < min_dist)
				{
					min_dist = dist;
					closest_tPos = tPos;
					bFound = true;
				}
			}

			if(bFound)
			{
				// The sum of the closest target positions
				sum_tPos += closest_tPos;
				//std::cout << sum_tPos << std::endl;
			}
			else
			{
				// Reduce the num because no correspondence found
				numSourceNeighbors--;
			}
		}

		// No correspondence found, so we set zero
		if(numSourceNeighbors == 0)
		{
			this->m_SurfaceVertexToCorrespondenceMapContainer->InsertElement(id,zeroPos);
			continue;
		}

		// Checks
		if(numSourceNeighbors < 0)
			itkExceptionMacro("Invalid number of source neighbor vertices : "<< numSourceNeighbors);

		if(sum_tPos == zeroPos)
			itkExceptionMacro("Invalid result vector!");

		// s'
		sum_tPos /= (double)numSourceNeighbors;

		// s
		sPos = this->m_FEMObject->GetNode(id)->GetCoordinates();

		// s'-s
		distPos = sum_tPos - sPos;

		// vertex normal
		normal = this->m_SurfaceVertexToNormalMapContainer->ElementAt(id);
		//std::cout << "normal : " << normal << std::endl;

		// project s'-s on the normal
		distPos_proj = dot_product(distPos,normal) * normal;

		// 1. Evaluate the correspondence based on the edges length
		EvaluateCorrespondenceWithEdgeLength(distPos_proj,id,this->GetScaleFactorForEdges());

		// Compute the projection of s'
		sum_tPos_proj = distPos_proj + sPos;

		// 2. Evaluate the correspondence based on the tets volumes
		EvaluateCorrespondenceWithVolume(distPos_proj,sum_tPos_proj,id,this->GetScaleFactorForVolumes());

		// Insert the projected vector displacement to the map
		this->m_SurfaceVertexToCorrespondenceMapContainer->InsertElement(id,distPos_proj);
	}
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::EvaluateCorrespondenceWithEdgeLength(VectorType& dVec,unsigned int vId,const double fac)
 {
	if(fac > 1.0 || fac <= 0.0)
		itkExceptionMacro("Invalid scale factor for edges:" << fac);

	typedef Element::Node NodeType;
	NodeType::Pointer pNode = this->m_FEMObject->GetNode(vId);
	if(pNode.IsNull())
		itkExceptionMacro("Node is NULL!");

	NodeType::SetOfElements adjacentElements = pNode->m_elements;
	if(adjacentElements.size() < 2)
		itkExceptionMacro("The number of adjacent tetrahedra is " << adjacentElements.size());

	NodeType::SetOfElements::iterator iter;
	Element* pTetr = NULL;
	double averageEdgeLength = 0;
	int numEdges = 0;

	for(iter = adjacentElements.begin(); iter != adjacentElements.end(); iter++)
	{
		pTetr = *iter;
		if(!pTetr)
			itkExceptionMacro("The tetr of the adjacent list is NULL!");

		std::vector<std::vector<int> > edgeIds = pTetr->GetEdgeIds();
		if(edgeIds.size() != 6)
			itkExceptionMacro("The tetrahedron doesn't have 6 edges!");

		unsigned int startId,endId;
		VectorType startPos,endPos;
		VectorType distVec;
		for(unsigned int i=0; i<6; i++)
		{
			startId = edgeIds.at(i).at(0);
			endId = edgeIds.at(i).at(1);

			startPos = pTetr->GetNode(startId)->GetCoordinates();
			endPos = pTetr->GetNode(endId)->GetCoordinates();

			distVec = endPos-startPos;

			averageEdgeLength += distVec.magnitude();
			numEdges++;
		}
	}
	if(numEdges < 1)
		itkExceptionMacro("Number of adjacent edges < 1");

	averageEdgeLength /= (double)numEdges;
	averageEdgeLength *= fac; // we scale the average length

	double dNorm = dVec.magnitude();

	// scaling
	double ratio = averageEdgeLength/dNorm;
	if(ratio < 1.0)
	{
		//std::cout << dVec << std::endl;
		dVec *= ratio;
		//std::cout << dVec << std::endl;
	}
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::EvaluateCorrespondenceWithVolume(VectorType& dVec,VectorType& newPos,unsigned int vId,const double fac)
 {
	const double ZERO6 = 0.000001;

	if(fac > 1.0 || fac <= 0.0)
		itkExceptionMacro("Invalid scale factor for Volume :" << fac);

	if(fabs(fac - 1.0) < 0.1)
		return;

	typedef Element::Node NodeType;
	NodeType::Pointer pNode = this->m_FEMObject->GetNode(vId);
	if(pNode.IsNull())
		itkExceptionMacro("Node is NULL!");

	NodeType::SetOfElements adjacentElements = pNode->m_elements;
	if(adjacentElements.size() < 3)
		itkExceptionMacro("The number of adjacent tetrahedra is " << adjacentElements.size());

	NodeType::SetOfElements::iterator iter;
	Element* pTetr = NULL;
	double volume = 0;
	std::vector<VectorType> tetrCoord;
	tetrCoord.resize(4);

	double dNorm = dVec.magnitude();
	bool bFoundNegativeVolume;
	const double minNorm = 0.001;

	int count = 0;
	iter = adjacentElements.begin();

	for(;;)
	{
		bFoundNegativeVolume = false;
		while(iter != adjacentElements.end())
		{
			pTetr = *iter;
			if(!pTetr)
				itkExceptionMacro("The tetr of the adjacent list is NULL!");

			// Found the coordinates of the new tet
			for(unsigned int i=0; i<4; i++)
			{
				if(pTetr->GetNode(i) == pNode)
					tetrCoord.at(i) = newPos;
				else
					tetrCoord.at(i) = pTetr->GetNode(i)->GetCoordinates();
			}

			volume = CalculateTetrahedronVolume(tetrCoord.at(0),tetrCoord.at(1),tetrCoord.at(2),tetrCoord.at(3));
			if(volume < ZERO6)
			{
				bFoundNegativeVolume = true;
				break;
			}
			iter++;
			count++;
		}

		// Reduce the displacement vector by a factor
		if(bFoundNegativeVolume)
		{
			dVec *= fac;
			dNorm = dVec.magnitude();

			// In case where the norm is zero we apply zero displacements and the new position
			// equals to the undeformed position
			if(dNorm < minNorm)
			{
				dVec.fill(0);
				newPos = pNode->GetCoordinates();
				return;
			}

			//std::cout << "newPos_before_reduction : " << newPos << std::endl;
			newPos = dVec + pNode->GetCoordinates();
			//std::cout << "newPos_after_reduction : " << newPos << std::endl;
			//			std::cout << "Initial Tetr : " << std::endl;
			//			std::cout << pTetr->GetNode(0)->GetCoordinates() << std::endl;
			//			std::cout << pTetr->GetNode(1)->GetCoordinates() << std::endl;
			//			std::cout << pTetr->GetNode(2)->GetCoordinates() << std::endl;
			//			std::cout << pTetr->GetNode(3)->GetCoordinates() << std::endl;
			//			std::cout << "New Tetr before reduction : " << std::endl;
			//			std::cout << tetrCoord.at(0) << std::endl;
			//			std::cout << tetrCoord.at(1) << std::endl;
			//			std::cout << tetrCoord.at(2) << std::endl;
			//			std::cout << tetrCoord.at(3) << std::endl;
		}
		else
			break;
	}
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::AssembleF()
 {
	//std::cout << "Assembling Forces..." << std::endl;
	const double pointTensorPonderation = GetLandmarkTensorPonderation();

	this->GetLinearSystemWrapper()->InitializeVector(m_ForceIndex);

	LoadContainerType *container = this->m_FEMObject->GetLoadContainer();
	if(!container)
		itkExceptionMacro("Missing container");

	LoadContainerIterator it = container->Begin();

	for(;it != container->End(); ++it)
	{
		Load::Pointer load = it.Value();
		LoadNoisyLandmark* landmark = dynamic_cast<LoadNoisyLandmark*>(load.GetPointer());
		if(!landmark)
			itkExceptionMacro("Landmark is NULL");

		if(landmark->IsOutlier())
			itkExceptionMacro("Found an Outlier Landmark!");

		const double confidence = landmark->GetConfidence();
		const VectorType & realDisplacement = landmark->GetRealDisplacement();
		const MatrixType & tensor = landmark->GetLandmarkTensor();
		const VectorType & shape = landmark->GetShape();
		const Element * element = landmark->GetContainedElement();

		const FEMIndexType numberOfDOFs = element->GetNumberOfDegreesOfFreedomPerNode();
		const FEMIndexType numberOfNodes = element->GetNumberOfNodes();
		float tradeOff = this->GetLambda();

		for(FEMIndexType m = 0;m < numberOfNodes; ++m)
		{
			double barCoor = shape[m];
			const VectorType weightedRealDisplacement = confidence * barCoor * pointTensorPonderation * tradeOff * ((tensor) * realDisplacement);

			for(FEMIndexType j = 0; j < numberOfDOFs; ++j)
			{
				const int index = element->GetDegreeOfFreedom(m * numberOfDOFs + j);
				this->GetLinearSystemWrapper()->AddVectorValue(index, weightedRealDisplacement[j], this->m_ForceIndex);
			}
		}
	}
}


template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::AssembleGlobalMatrixFromLandmarksAndMeshMatrices()
 {
	//std::cout << "Assembling Global Matrix From Landmarks And Mesh Matrices..." << std::endl;
	this->m_ls->CopyMatrix( this->m_MeshStiffnessMatrixIndex, this->m_StiffnessMatrixIndex );
	this->m_ls->AddMatrixMatrix( this->m_StiffnessMatrixIndex, this->m_LandmarkStiffnessMatrixIndex);
}

template <class TPixel,unsigned int VDimension>
float CompressorSolver<TPixel,VDimension>
::GetLandmarkTensorPonderation(void) const
{
  const LoadContainerType * loadContainer = this->m_FEMObject->GetLoadContainer();

  if(!loadContainer)
    {
    itkExceptionMacro("Missing load container");
    }

  const NodeContainerType * nodeContainer = this->m_FEMObject->GetNodeContainer();

  if(!nodeContainer)
    {
    itkExceptionMacro("Missing node container");
    }

  const float ponderation =
    static_cast<float>( nodeContainer->Size() ) /
    static_cast<float>( loadContainer->Size() );

  return ponderation;
}


template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::AssembleLandmarkStiffnessMatrix()
 {
	//std::cout << "Assembling Landmarks Stiffness Matrix..." << std::endl;

	// Assemble the contribution matrix of the landmarks
	const double pointTensorPonderation = this->GetLandmarkTensorPonderation();
	LoadContainerType *container = this->m_FEMObject->GetLoadContainer();
	if(!container)
	{
		itkExceptionMacro("Missing container");
	}

	LoadContainerIterator it = container->Begin();

	for(;it != container->End(); ++it)
	{
		Load::Pointer load = it.Value();

		LoadNoisyLandmark *landmark = dynamic_cast<LoadNoisyLandmark*>(load.GetPointer());
		if(landmark == NULL)
			itkExceptionMacro("Encounter landmark that is not a LoadNoisyLandmark");

		if(landmark->IsOutlier())
			itkExceptionMacro("Found Outlier Landmark!");

		const double confidence = landmark->GetConfidence();
		const MatrixType & tens = landmark->GetLandmarkTensor();
		const VectorType & shape = landmark->GetShape();
		const Element * element = landmark->GetContainedElement();
		const FEMIndexType numberOfDOFs = element->GetNumberOfDegreesOfFreedomPerNode();
		const FEMIndexType numberOfNodes = element->GetNumberOfNodes();
		float tradeOff = this->GetLambda();

		// fill the diagonal matrices
		for(FEMIndexType k = 0; k < numberOfNodes; ++k)
		{
			const double barCoor = shape[k] * shape[k];

			for(FEMIndexType n = 0; n < numberOfDOFs; n++)
			{
				for(FEMIndexType m = 0; m < numberOfDOFs; m++)
				{
					const int dofn = element->GetDegreeOfFreedom(k * numberOfDOFs + n);
					const int dofm = element->GetDegreeOfFreedom(k * numberOfDOFs + m);
					const float value = static_cast<float>( barCoor * tradeOff * pointTensorPonderation * (tens(n,m)) * confidence );

					this->m_ls->AddMatrixValue(dofn, dofm, value, m_LandmarkStiffnessMatrixIndex);
				}
			}
		}

		// fill the extradiagonal matrices
		for(FEMIndexType i = 0; i < numberOfNodes - 1; i++)
		{
			for(FEMIndexType j = i + 1; j < numberOfNodes; j++)
			{
				const double barCoor = shape[i] * shape[j];

				for(FEMIndexType n = 0; n < numberOfDOFs; n++)
				{
					for(FEMIndexType m = 0; m < numberOfDOFs; m++)
					{
						const int dofn = element->GetDegreeOfFreedom(i * numberOfDOFs + n);
						const int dofm = element->GetDegreeOfFreedom(j * numberOfDOFs + m);
						const float value = static_cast<float>( barCoor * tradeOff * pointTensorPonderation * (tens(n, m)) * confidence );

						this->m_ls->AddMatrixValue(dofn, dofm, value, m_LandmarkStiffnessMatrixIndex);
						this->m_ls->AddMatrixValue(dofm, dofn, value, m_LandmarkStiffnessMatrixIndex);

					}
				}
			}
		}
	}
}


template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::ComputeLandmarksTensor()
 {
	//std::cout << "Computing Landmarks Tensor..." << std::endl;

	// Compute mesh Surface Vertices tensor weighted by a structure tensor if exists
	LoadContainerType *container = this->m_FEMObject->GetLoadContainer();
	if(!container)
		itkExceptionMacro("Missing container");

	LoadContainerIterator it = container->Begin();

	for(; it != container->End(); ++it)
	{
		Load::Pointer pLoad = it.Value();

		if(pLoad.IsNull())
			itkExceptionMacro("Load is NULL!");
		LoadNoisyLandmark *pLandmark = dynamic_cast<LoadNoisyLandmark*>(pLoad.GetPointer());
		if(!pLandmark)
			itkExceptionMacro("Landmark is NULL!");

		// We do not accept outliers
		if(pLandmark->IsOutlier())
			itkExceptionMacro("Found outlier landmark!");

		const VectorType & shape = pLandmark->GetShape();
		const Element * element = pLandmark->GetContainedElement();
		if(!element)
			itkExceptionMacro("Element is NULL!");

		const FEMIndexType numberOfDOFs = element->GetNumberOfDegreesOfFreedomPerNode();
		const FEMIndexType numberOfNodes = element->GetNumberOfNodes();

		MatrixType nodeTensor(numberOfDOFs, numberOfDOFs);
		MatrixType landmarkTensor(numberOfDOFs, numberOfDOFs);

		landmarkTensor.fill(0.0);

		for(FEMIndexType nodeId = 0; nodeId < numberOfNodes; ++nodeId)
		{
			for(FEMIndexType dofXId = 0; dofXId < numberOfDOFs; dofXId++)
			{
				for(FEMIndexType dofYId = 0; dofYId < numberOfDOFs; dofYId++)
				{
					unsigned nx = element->GetDegreeOfFreedom(nodeId * numberOfDOFs + dofXId);
					unsigned ny = element->GetDegreeOfFreedom(nodeId * numberOfDOFs + dofYId);
					nodeTensor[dofXId][dofYId] = this->m_ls->GetMatrixValue(nx, ny, m_MeshStiffnessMatrixIndex);
				}
			}
			landmarkTensor += nodeTensor * shape[nodeId];
		}
		if(pLandmark->HasStructureTensor())
			landmarkTensor *= pLandmark->GetStructureTensor();

		pLandmark->SetLandmarkTensor(landmarkTensor);
	}
}


template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::AssembleMeshStiffnessMatrix()
{
	//std::cout << "Assembling Mesh Stiffness Matrix..." << std::endl;
 // Assemble the mechanical stiffness matrix from the mesh

  // if no DOFs exist in a system, we have nothing to do
  if( this->m_NGFN <= 0 )
    {
    return;
    }

  // assemble the mechanical matrix by stepping over all elements
  FEMIndexType numberOfElements = this->m_FEMObject->GetNumberOfElements();
  for( FEMIndexType i = 0; i < numberOfElements; i++)
    {
    // call the function that actually moves the element matrix to the master matrix.
    Element::Pointer element = this->m_FEMObject->GetElement(i);

    this->AssembleElementMatrixWithID(element, this->m_MeshStiffnessMatrixIndex);
    }
}

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::AssembleElementMatrixWithID(const Element::Pointer & element, unsigned int matrixIndex)
{
  // copy the element stiffness matrix for faster access.
  Element::MatrixType Ke;

  element->GetStiffnessMatrix(Ke);

  // same for number of DOF
  const FEMIndexType numberOfDOFs = element->GetNumberOfDegreesOfFreedom();

  // step over all rows in element matrix
  for( FEMIndexType j = 0; j < numberOfDOFs; j++ )
    {
    // step over all columns in element matrix
    for( FEMIndexType k = 0; k < numberOfDOFs; k++ )
      {
      // error checking. all GFN should be >= 0 and < NGFN
      const unsigned int dofj = element->GetDegreeOfFreedom(j);
      const unsigned int dofk = element->GetDegreeOfFreedom(k);
      if( dofj >= this->m_NGFN || dofk >= this->m_NGFN )
        {
        throw FEMExceptionSolution(__FILE__, __LINE__, "Solver::AssembleElementMatrix()", "Illegal GFN!");
        }

      if( Ke[j][k] != Float(0.0) )
        {
        this->GetLinearSystemWrapper()->AddMatrixValue(dofj, dofk, Ke[j][k], matrixIndex);
        }
      }
    }
}

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::Initialization()
 {
	this->SetLinearSystemWrapper(&m_Itpack);
	const FEMIndexType maximumNonZeroMatrixEntriesFactor = 100;
	const FEMIndexType maxNumberOfNonZeroValues = this->m_NGFN * maximumNonZeroMatrixEntriesFactor;
	if( maxNumberOfNonZeroValues > NumericTraits< FEMIndexType >::max() / 2 )
	{
		itkExceptionMacro("Too large system of equations");
	}

	this->m_Itpack.SetMaximumNonZeroValuesInMatrix( maxNumberOfNonZeroValues );

	// the NGFN is determined once the FEMObject is finalized
	this->GetLinearSystemWrapper()->SetSystemOrder(this->m_NGFN);
	this->GetLinearSystemWrapper()->SetNumberOfVectors(1);
	this->GetLinearSystemWrapper()->SetNumberOfSolutions(1);
	this->GetLinearSystemWrapper()->SetNumberOfMatrices(3);
	this->GetLinearSystemWrapper()->InitializeMatrix(m_MeshStiffnessMatrixIndex);
	this->GetLinearSystemWrapper()->InitializeMatrix(m_LandmarkStiffnessMatrixIndex);
	this->GetLinearSystemWrapper()->InitializeMatrix(m_StiffnessMatrixIndex);
	this->GetLinearSystemWrapper()->InitializeVector(m_ForceIndex);
	this->GetLinearSystemWrapper()->InitializeSolution(m_SolutionIndex);
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::AddSurfaceVertexToItsOwnNeighborsSet()
 {
	if(!this->m_SurfaceVertexToLabelMapContainer)
		itkExceptionMacro("Surface Vertex To Label Map Container is NULL!");
	if(!this->m_SurfaceVertexToNeighborVerticesMapContainer)
		itkExceptionMacro("Surface Vertex To Neighbor Vertices Map Container is NULL!");

	typename SurfaceVertexToLabelMapContainerType::Iterator vertexToLabelIt;
	unsigned int id = 0;
	// Loop over the surface vertices
	for(vertexToLabelIt = this->m_SurfaceVertexToLabelMapContainer->Begin(); vertexToLabelIt != this->m_SurfaceVertexToLabelMapContainer->End();  vertexToLabelIt++)
	{
		id = vertexToLabelIt->Index();
		// Get the neighbors of the vertex
		std::set<unsigned int>* pNeighbors = this->m_SurfaceVertexToNeighborVerticesMapContainer->ElementAt(id);
		if(!pNeighbors)
			itkExceptionMacro("The neighbors vertices vector is NULL!");
		if(pNeighbors->size() < 3)
			itkExceptionMacro("The size of the neighbor set is less than 3!");

		// Add the id of the current vertex
		//std::cout << pNeighbors->size() << std::endl;
		pNeighbors->insert(id);
		//std::cout << pNeighbors->size() << std::endl;
	}
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::RemoveSurfaceVertexFromItsOwnNeighborsSet()
 {
	if(!this->m_SurfaceVertexToLabelMapContainer)
		itkExceptionMacro("Surface Vertex To Label Map Container is NULL!");
	if(!this->m_SurfaceVertexToNeighborVerticesMapContainer)
		itkExceptionMacro("Surface Vertex To Neighbor Vertices Map Container is NULL!");

	typename SurfaceVertexToLabelMapContainerType::Iterator vertexToLabelIt;
	unsigned int id = 0;
	// Loop over the surface vertices
	for(vertexToLabelIt = this->m_SurfaceVertexToLabelMapContainer->Begin(); vertexToLabelIt != this->m_SurfaceVertexToLabelMapContainer->End();  vertexToLabelIt++)
	{
		id = vertexToLabelIt->Index();
		// Get the neighbors of the vertex
		std::set<unsigned int>* pNeighbors = this->m_SurfaceVertexToNeighborVerticesMapContainer->ElementAt(id);
		if(!pNeighbors)
			itkExceptionMacro("The neighbors vertices vector is NULL!");
		if(pNeighbors->size() < 3)
			itkExceptionMacro("The size of the neighbor set is less than 3!");

		// Remove the id of the current vertex
		pNeighbors->erase(id);
	}
 }

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::PrintForces(const char* fName)
 {
	if(fName == NULL)
		return;

	FILE* pFile = fopen(fName,"w");
	if(pFile == NULL)
		return;

	std::cout <<  "Printing Forces at : " << fName << std::endl;

	Float val = 0;
	int NumDofs = this->m_ls->GetSystemOrder();
	std::cout <<  "Num DOF's : " << NumDofs << std::endl;
	for( int i = 0; i < NumDofs; i++ )
	{
		val = this->m_ls->GetVectorValue(i,m_ForceIndex);
		fprintf(pFile,"%6.8f\n",val);
	}
	fclose(pFile);
	std::cout << "Done" << std::endl;
}

template <class TPixel,unsigned int VDimension>
void CompressorSolver<TPixel,VDimension>
::PrintSolution(const char * fName)
{
	if(fName == NULL)
		return;

	FILE* pFile = fopen(fName,"w");
	if(pFile == NULL)
		return;

	std::cout <<  "Printing Solution Displacements at : " << fName << std::endl;

	Float val = 0;
	int NumDofs = this->m_ls->GetSystemOrder();
	std::cout <<  "Num DOF's : " << NumDofs << std::endl;
	for( int i = 0; i < NumDofs; i++ )
	{
		// Get the solution displacements
		val = this->m_ls->GetSolutionValue(i,m_SolutionIndex);
		fprintf(pFile,"%10.10f\n",val);
	}
	fclose(pFile);
	std::cout << "Done" << std::endl;
}

}
}

#endif /* ITKFEMCOMPRESSORSOLVER_HXX_ */
