/*
 * itkFEMCompressor.hxx
 *
 */

#ifndef ITKFEMCOMPRESSOR_HXX_
#define ITKFEMCOMPRESSOR_HXX_

#include "itkImageDuplicator.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkFEMCompressorSolver.h"
#include "vtkNew.h"

namespace itk
{
namespace fem
{

template<class TInputMesh, class TInputImage, class TOutputMesh>
FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::FEMCompressor()
 {
	this->m_NonConnectivityType = VERTEX_NON_CONNECTIVITY;
	this->m_Lambda = 1.0;
	this->m_Image = NULL;
	this->m_SurfaceVertexToLabelMapContainer = NULL;
	this->m_SurfaceVertexToNeighborVerticesMapContainer = NULL;
	this->m_SurfaceTrianglesContainer = NULL;
	this->m_TargetPointSet = NULL;
	this->m_FEMObject = NULL;
	this->m_FEMCompressorSolver = NULL;
	this->m_NumIterations = 5;
	this->m_MaterialBrain = NULL;
	this->m_YoungModulus = 2100.0 * 1e-06; // YoungModulus (2100 Pa = 2100 * 10^{-6} N/mm2)
	this->m_PoissonRatio = 0.45;
	this->m_TetrahedraDataContainer = NULL;
	this->m_SelectionMapImage = NULL;
	this->m_OutputSourcePointsFileName = "";
	this->m_OutputTargetPointsFileName = "";
	this->m_OutputSurfaceTrianglesFileName = "";
 }


template<class TInputMesh, class TInputImage, class TOutputMesh>
FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::~FEMCompressor()
{

}


template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::GenerateData()
 {
	if(!this->GetInput())
		itkExceptionMacro("The input mesh is NULL." << "SetInput() must be called first.");

	// Check the input mesh for inverted elements
	if(!CheckTetrahedraVolumes(this->GetInput()))
		itkExceptionMacro("The input Mesh Contains inverted elements!");

	// Extraction of the surface vertices ids
	ExtractSourcePoints();

	// Compute the neighbors of surface vertices
	BuildSourcePointsConnectivity();

	// The target points are extracted from the label image and stored in a point set.
	ExtractTargetPoints();

	// Create and initialize the FEM Object
	InitializeFEMObject();

	// Set parameters and Update the Compressor solver
	this->m_FEMCompressorSolver = FEMCompressorSolverType::New();
	this->m_FEMCompressorSolver->SetInput( this->m_FEMObject );
	this->m_FEMCompressorSolver->SetSurfaceVertexToLabelMapContainer(this->m_SurfaceVertexToLabelMapContainer);
	this->m_FEMCompressorSolver->SetSurfaceVertexToNeighborVerticesMapContainer(this->m_SurfaceVertexToNeighborVerticesMapContainer);
	this->m_FEMCompressorSolver->SetSurfaceTrianglesContainer(this->m_SurfaceTrianglesContainer);
	this->m_FEMCompressorSolver->SetTargetPointSet(this->m_TargetPointSet);
	this->m_FEMCompressorSolver->SetLambda(this->GetLambda());
	this->m_FEMCompressorSolver->SetNumIterations(this->GetNumIterations());
	this->m_FEMCompressorSolver->SetTetrahedraDataContainer(this->m_TetrahedraDataContainer);// for printing tets with labels
	this->m_FEMCompressorSolver->Update();

	UpdateOutputMesh();
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::UpdateOutputMesh()
 {
	this->CopyInputMeshToOutputMeshPoints();
	this->CopyInputMeshToOutputMeshPointData();
	this->CopyInputMeshToOutputMeshCellLinks();
	this->CopyInputMeshToOutputMeshCells();
	this->CopyInputMeshToOutputMeshCellData();

	// Copy the FEMObject nodes coordinates to the output mesh
	FEMObjectPointerType pFEMObject = this->m_FEMCompressorSolver->GetInput();
	//	FEMObjectPointerType pFEMObject = this->m_FEMCompressorSolver->GetOutput();
	NodeContainerType* pNodeContainer =  pFEMObject->GetNodeContainer();
	MeshType* pOutputMesh = this->GetOutput();
	if(pOutputMesh->GetNumberOfPoints() != pNodeContainer->Size())
		itkExceptionMacro("Output mesh and output FEMObject have different number of points!");

	typename NodeContainerType::Iterator iter;
	VectorType pos;
	PointType pt;
	for(iter = pNodeContainer->Begin(); iter != pNodeContainer->End();  iter++)
	{
		pos.SetVnlVector(iter.Value()->GetCoordinates());
		for(int i=0; i<MeshDimension; i++)
			pt[i] = pos[i];

		pOutputMesh->SetPoint(iter.Index(),pt);
	}
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::ExtractTargetPoints()
 {
	std::cout << "" << std::endl;
	std::cout << "Extracting  target Points..." << std::endl;

	ImageType* pInputImage = this->GetImage();
	if(!pInputImage)
		itkExceptionMacro("The input Image is NULL!");

	this->ComputeConnectivityOffsets();

	// Create the selection map Image
	this->m_SelectionMapImage = SelectionMapImageType::New();
	ImageRegionType inputRegion = pInputImage->GetLargestPossibleRegion();
	this->m_SelectionMapImage->SetRegions(inputRegion);
	this->m_SelectionMapImage->SetOrigin(pInputImage->GetOrigin());
	this->m_SelectionMapImage->SetSpacing(pInputImage->GetSpacing());
	this->m_SelectionMapImage->SetDirection(pInputImage->GetDirection());
	this->m_SelectionMapImage->Allocate();
	const typename SelectionMapImageType::PixelType zero = 0;
	this->m_SelectionMapImage->FillBuffer(zero);
	this->m_SelectionMapImage->Update();

	// copy labeled image into selectionMap
	itk::ImageRegionConstIterator< ImageType> maskItr( pInputImage, inputRegion);
	itk::ImageRegionIterator< SelectionMapImageType > mapItr( this->m_SelectionMapImage, inputRegion );
	const typename SelectionMapImageType::PixelType one = 1;
	for ( maskItr.GoToBegin(), mapItr.GoToBegin(); !maskItr.IsAtEnd(); ++maskItr, ++mapItr )
	{
		if(maskItr.Get() != 0)
			mapItr.Set(one);
	}

	// Collect the label of the tissues
	typedef std::map<ImagePixelType,unsigned int> LabelMapType;
	LabelMapType labelMap;

	typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
	ConstIteratorType inputIt(pInputImage, pInputImage->GetLargestPossibleRegion());
	ImagePixelType label = 0;
	while( !inputIt.IsAtEnd() )
	{
		label = inputIt.Get();
		if(labelMap.find(label) != labelMap.end())
			labelMap[label]++;
		else
			labelMap[label] = 1;

		++inputIt;
	}

	typename LabelMapType::iterator it;
	unsigned int numVoxels = 0;
	for(it = labelMap.begin(); it != labelMap.end(); it++)
	{
		numVoxels += it->second;
		std::cout << "Label " << it->first << " has " << it->second << " voxels" << std::endl;
	}
	if(numVoxels != pInputImage->GetLargestPossibleRegion().GetSize()[0]*
			pInputImage->GetLargestPossibleRegion().GetSize()[1]*
			pInputImage->GetLargestPossibleRegion().GetSize()[2])
		itkExceptionMacro("Invalid number of total image voxels!");
	std::cout << "Number of input voxels : " << numVoxels << std::endl;

	if(labelMap.size() != 2)
        itkExceptionMacro("Invalid input image! " << labelMap.size() << " labels detected (including the background). "
                "Currently the method supports only 2 labels (background and tissue).");

	// Voxels in neighborhood of the current voxel
	ImagePixelType neighborLabel;
	ImageIndexType neighborIndex;
	ImageIndexType index;

	unsigned int numTargetVoxels = 0;

	typedef std::map<PixelTypeSet,std::vector<PointType>* > Label2MapType;
	Label2MapType label2Map;
	typename Label2MapType::iterator label2MapIter;
	std::vector<PointType>* pPointsVector = NULL;
	PointType pt;

	std::vector<PixelTypeSet> labelSetVector;

	// Loop of the labeled tissues
	for(it = labelMap.begin(); it != labelMap.end(); it++)
	{
		// The current label
		label = it->first;

		// Skip zero labels (background)
		if(label == 0)
			continue;

		inputIt.GoToBegin();

		// Loop over the image
		while( !inputIt.IsAtEnd() )
		{
			// The voxel label should be equal to the curent tissue label
			if(inputIt.Get() != label)
			{
				++inputIt;
				continue;
			}

			index = inputIt.GetIndex();

			// The index should be inside the selection map (PixelVal = 1)
			if(this->m_SelectionMapImage->GetPixel(index) == zero)
			{
				++inputIt;
				continue;
			}

			// Insert the current label
			PixelTypeSet labelSet;
			labelSet.insert(label);

			// We check in x,y,z directions with offset +-1 voxels from the current voxel.
			// We do not need to check the center voxel.
			for(unsigned int Dim = 0; Dim < ImageSizeType::Dimension; Dim++)
			{
				neighborIndex = index;

				for(int i = index[Dim]-1; i < index[Dim]+2; i = i+2)
				{
					neighborIndex[Dim] = i;

					// Test if the index is inside the image
					if(!inputRegion.IsInside(neighborIndex))
						continue;

					// The neighbor label
					neighborLabel = pInputImage->GetPixel(neighborIndex);
					labelSet.insert(neighborLabel);
				}
			}

			// If the edge point is boundary label set size will be larger than one
			if(labelSet.size() > 1)
			{
				// Check if the label set exists
				bool bSame = false;
				for(unsigned int j=0; j<labelSetVector.size(); j++)
				{
					if(SameLabelSets(labelSetVector.at(j),labelSet))
					{
						bSame = true;
						break;
					}
				}

				if(!bSame)
				{
					// Convert the center index to point coordinates
					pInputImage->TransformIndexToPhysicalPoint(index,pt);

					label2MapIter = label2Map.find(labelSet);

					// The label set exist
					if(label2MapIter != label2Map.end())
					{
						pPointsVector = (*label2MapIter).second;
					}
					// The label set doesn't exist
					else
					{
						pPointsVector = new std::vector<PointType>();
						label2Map[labelSet] = pPointsVector;
					}

					pPointsVector->push_back(pt);
					numTargetVoxels++;
					if(this->m_NonConnectivityType != NO_NON_CONNECTIVITY)
						MarkOffConnectedPoints(index);
				}
			}
			++inputIt;
		}

		// Add the label sets to the check vector
		labelSetVector.clear();
		for(label2MapIter = label2Map.begin(); label2MapIter != label2Map.end(); label2MapIter++)
			labelSetVector.push_back(label2MapIter->first);
	}


	// Create the target point set
	this->m_TargetPointSet = TargetPointSetType::New();
	typename TargetPointSetType::PointIdentifier id = 0;
	for(label2MapIter = label2Map.begin(); label2MapIter != label2Map.end(); label2MapIter++)
	{
		PixelTypeSet labelSet = label2MapIter->first;
		std::vector<PointType>* pPoints = label2MapIter->second;
		if(!pPoints)
			itkExceptionMacro("Invalid vector of points!");

		for(unsigned int i = 0; i < pPoints->size(); i++)
		{
			this->m_TargetPointSet->SetPoint(id, pPoints->at(i));
			this->m_TargetPointSet->SetPointData(id, labelSet);
			id++;
		}
	}

	// Check
	if(this->m_TargetPointSet->GetPoints()->Size() == 0)
		itkExceptionMacro("No target points were extracted from the labeled image!");

	// Prints
	std::cout << "" << std::endl;
	for(label2MapIter = label2Map.begin(); label2MapIter != label2Map.end(); label2MapIter++)
	{
		typename PixelTypeSet::iterator itt;
		PixelTypeSet labelSet = label2MapIter->first;
		std::cout << "label set:";
		for(itt = labelSet.begin(); itt != labelSet.end(); ++itt)
			std::cout << " " << *itt;
		std::cout << std::endl;
		std::cout << " has " << label2MapIter->second->size() << " voxels " << std::endl;
	}

	std::cout << "Number of target Points : " << numTargetVoxels << std::endl;

	// Prints
	if(!this->GetOutputTargetPointsFileName().empty())
		PrintTargetPointsVTK(this->m_TargetPointSet,this->GetOutputTargetPointsFileName());

	//PrintTargetPointsITK(this->m_TargetPointSet,pInputImage,"targetPoints.nii.gz");
 }


template<class TInputMesh, class TInputImage, class TOutputMesh>
void
FEMCompressor< TInputMesh, TInputImage, TOutputMesh >
::MarkOffConnectedPoints(ImageIndexType index)
{
	// Mark off connected points
	const typename SelectionMapImageType::PixelType ineligeblePoint = 0;
	for ( SizeValueType j = 0, n = this->m_NonConnectivityOffsets.size(); j < n; j++ )
	{
		ImageIndexType idx = index;
		idx += m_NonConnectivityOffsets[ j ];
		this->m_SelectionMapImage->SetPixel( idx, ineligeblePoint );
	}
}

template<class TInputMesh, class TInputImage, class TOutputMesh>
void
FEMCompressor< TInputMesh, TInputImage, TOutputMesh >
::ComputeConnectivityOffsets( )
 {
	if ( this->m_NonConnectivityType == VERTEX_NON_CONNECTIVITY ||
		 this->m_NonConnectivityType == EDGE_NON_CONNECTIVITY ||
		 this->m_NonConnectivityType == FACE_NON_CONNECTIVITY )
	{
		if(this->m_NonConnectivityType  == VERTEX_NON_CONNECTIVITY)
			std::cout << "Non-Connectivity Type : VERTEX_NON_CONNECTIVITY" << std::endl;
		else if(this->m_NonConnectivityType == EDGE_NON_CONNECTIVITY)
			std::cout << "Non-Connectivity Type : EDGE_NON_CONNECTIVITY" << std::endl;
		else
			std::cout << "Non-Connectivity Type : FACE_NON_CONNECTIVITY" << std::endl;

		this->m_NonConnectivityOffsets.clear();
		// use Neighbourhood to compute all offsets in radius 1
		Neighborhood< unsigned int, MeshDimension > neighborhood;
		neighborhood.SetRadius( NumericTraits< SizeValueType >::One );
		for ( SizeValueType i = 0, n = neighborhood.Size(); i < n; i++ )
		{
			OffsetType off = neighborhood.GetOffset( i );

			// count 0s offsets in each dimension
			unsigned numberOfZeros = 0;
			for ( unsigned j = 0; j < MeshDimension; j++ )
			{
				if ( off[j] == 0 )
					numberOfZeros++;
			}

			if ( this->m_NonConnectivityType <= numberOfZeros && numberOfZeros < MeshDimension )
				this->m_NonConnectivityOffsets.push_back(off);
		}
		return;
	}
	else if(this->m_NonConnectivityType == NO_NON_CONNECTIVITY)
	{
		std::cout << "Non-Connectivity Type : NO_NON_CONNECTIVITY" << std::endl;
		return;
	}
	else
	{
		itkExceptionMacro( "Cannot use non-connectivity of value " << m_NonConnectivityType
				<< ", expected a value in the range 0.." << MeshDimension << "." );
	}
}


template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::ExtractSourcePoints()
 {
	// Check mesh
	const MeshType* pMesh = this->GetInput();
	if(!pMesh)
		itkExceptionMacro("The input Mesh is NULL.");
	if(!pMesh->GetCells())
		itkExceptionMacro("Cells container is NULL.");
	if(!pMesh->GetCellData())
		itkExceptionMacro("Cells Data container is NULL.");
	if(pMesh->GetCells()->Size() != pMesh->GetCellData()->Size())
		itkExceptionMacro("Num of Cells is different than the number of Cells Data.");

	// Build the cell links
	pMesh->BuildCellLinks();

	// The surface triangles container
	this->m_SurfaceTrianglesContainer = SurfaceTrianglesContainerType::New();
	this->m_SurfaceTrianglesContainer->Initialize();

	// The surface vertices container
	this->m_SurfaceVertexToLabelMapContainer = SurfaceVertexToLabelMapContainerType::New();
	this->m_SurfaceVertexToLabelMapContainer->Initialize();

	//Mesh cell iterator
	CellIteratorType cellIt = pMesh->GetCells()->Begin();
	CellIteratorType cellEnd = pMesh->GetCells()->End();
	CellDataIteratorType cellDataIt = pMesh->GetCellData()->Begin();
	CellDataIteratorType cellDataItEnd = pMesh->GetCellData()->End();

	ImagePixelType tetr_label=0,linked_tetr_label=0;
	TetrahedronType* pTetr = NULL;
	enum TRI_TYPE {EXTERIOR_SURFACE = 0, INTERIOR_SURFACE = 1, INTERIOR = 2};
	TRI_TYPE triType;
	unsigned int numExtSurfTriangles = 0;

	// Loop over all cells
	while( cellIt != cellEnd )
	{
		pTetr = (TetrahedronType*)cellIt.Value();
		if(!pTetr)
			itkExceptionMacro("Found NULL tetrahedron!");

		if(pTetr->GetType() != CellType::TETRAHEDRON_CELL)
			itkExceptionMacro("Invalide Cell Type! Currently only Tetrahedral cells are supported. Current element type is : " << pTetr->GetType());

		unsigned int tet_id = cellIt.Index();
		unsigned int linked_tet_id = 0;


		tetr_label = cellDataIt.Value();
		if(tetr_label <= 0)
			itkExceptionMacro("Found Tetr with invalid label!");


		if(pTetr->GetNumberOfFaces() != 4)
			itkExceptionMacro("A tetrahedron has not 4 faces!");

		typename TetrahedronType::FaceAutoPointer pFace;
		typename TetrahedronType::VertexAutoPointer pVertex;
		typename MeshType::CellAutoPointer pLinkedTetr;

		// Get faces
		for(unsigned int i=0; i < pTetr->GetNumberOfFaces(); i++)
		{
			pTetr->GetFace(i,pFace);
			if(pFace->GetNumberOfVertices() != 3)
				itkExceptionMacro("A triangle has not 3 vertices!");

			unsigned int vTriId[3] = {0,0,0};
			unsigned int vTetId[4] = {0,0,0,0};

			// Get the vertices of the face
			PointIdIterator pIdIter = ((TriangleType*)pFace.GetPointer())->PointIdsBegin();
			PointIdIterator pIdEnd = ((TriangleType*)pFace.GetPointer())->PointIdsEnd();
			unsigned int index=0;
			while( pIdIter != pIdEnd )
			{
				vTriId[index] = *pIdIter;
				pIdIter++;
				index++;
			}

			triType = EXTERIOR_SURFACE;

			// For each vertex
			for(unsigned int j=0; j<3; j++)
			{
				// Get the linked tets of the vertex
				std::set<unsigned long> vertex2Tets = pMesh->GetCellLinks()->ElementAt(vTriId[j]);
				std::set<unsigned long>::iterator it;
				if(vertex2Tets.size() < 2) // might be one??
					itkExceptionMacro("Invalid number of linked tets!");

				// Loop over the linked tets
				for(it = vertex2Tets.begin(); it != vertex2Tets.end(); it++)
				{
					linked_tet_id = *it;
					if(linked_tet_id == tet_id)// skip the current tet
						continue;

					if(!pMesh->GetCell(linked_tet_id,pLinkedTetr))
						itkExceptionMacro("The linked cell doesn't exist!");

					linked_tetr_label = pMesh->GetCellData()->GetElement(linked_tet_id);
					if(linked_tetr_label <= 0)
						itkExceptionMacro("Found linked Tetr with invalid label!");

					// Check if the linked tetrahedron contains the vertices of the face
					pIdIter = ((TetrahedronType*)pLinkedTetr.GetPointer())->PointIdsBegin();
					pIdEnd = ((TetrahedronType*)pLinkedTetr.GetPointer())->PointIdsEnd();
					index=0;
					while( pIdIter != pIdEnd )
					{
						vTetId[index] = *pIdIter;
						pIdIter++;
						index++;
					}

					// Check if the tetrahedron contains the vertices of triangle
					if(GetNumOfCommonIds(vTriId,vTetId) == 3)
					{
						// If the label of the linked tet is different, the triangle is SURFACE_INTERIOR
						if(tetr_label != linked_tetr_label)
							triType = INTERIOR_SURFACE;
						// If the label of the linked tet is the same, the triangle is INTERIOR
						else
							triType = INTERIOR;

						break;
					}
				}
				if(triType != EXTERIOR_SURFACE)
					break;
			}

			// All vertices are interior surface
			if(triType == INTERIOR_SURFACE)
			{
				itkExceptionMacro("An interface triangle was found. Currently only a single-tissue input mesh is supported.");
			}

			if(triType == EXTERIOR_SURFACE)
			{
				// Check the orientation of the surface triangle.
				// This function might change the order of the vertices ids of the surface triangle.
				CheckSurfaceTriangleOrientation(triType,vTriId,pTetr,(TetrahedronType*)pLinkedTetr.GetPointer(),tetr_label,linked_tetr_label);

				std::pair<typename SurfaceVertexToLabelMapContainerType::ElementIdentifier,typename SurfaceVertexToLabelMapContainerType::Element> vertexToLabels_pair;
				vertexToLabels_pair.second.insert(tetr_label);
				vertexToLabels_pair.second.insert(0);

				// For each vertex of the surface triangle
				VectorTriIdType verticesTriId;
				for(unsigned int k=0; k< MeshDimension; k++)
				{
					// Update the surface vertex id to label map
					vertexToLabels_pair.first = vTriId[k];

					// The index exist, so we update the set
					if(this->m_SurfaceVertexToLabelMapContainer->IndexExists(vTriId[k]))
					{
						PixelTypeSet pixelSet = this->m_SurfaceVertexToLabelMapContainer->GetElement(vTriId[k]);
						pixelSet.insert(linked_tetr_label);
						this->m_SurfaceVertexToLabelMapContainer->SetElement(vTriId[k],pixelSet);
					}
					// The index doesn't exist
					else
					{
						this->m_SurfaceVertexToLabelMapContainer->CastToSTLContainer().insert(vertexToLabels_pair);
					}
					verticesTriId[k] = vTriId[k];
				}

				// Update the surface triangles container
				this->m_SurfaceTrianglesContainer->CastToSTLContainer().push_back(verticesTriId);
			}

			if(triType == EXTERIOR_SURFACE)
				numExtSurfTriangles++;
		}
		cellIt++;
		cellDataIt++;
	}

	std::cout << "Num of Surface Triangles : " << this->m_SurfaceTrianglesContainer->Size() << std::endl;
	std::cout << "Num of Surface Triangles : " << numExtSurfTriangles << std::endl;
	std::cout << "Num of Surface Vertices : " << this->m_SurfaceVertexToLabelMapContainer->Size() << std::endl;
	std::cout << "Num of Source Points " << this->m_SurfaceVertexToLabelMapContainer->Size() << std::endl;

	// Prints
	if(!this->GetOutputSourcePointsFileName().empty())
		PrintSourcePointsVTK(this->m_SurfaceVertexToLabelMapContainer,this->GetInput()->GetPoints(),this->GetOutputSourcePointsFileName());
	if(!this->GetOutputSurfaceTrianglesFileName().empty())
		PrintSurfaceTriangles(this->m_SurfaceTrianglesContainer,this->GetInput()->GetPoints(),this->GetOutputSurfaceTrianglesFileName());
 }


template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::CheckSurfaceTriangleOrientation(int triType,unsigned int ids[3],TetrahedronType* pTetr,TetrahedronType* pLinkedTetr,ImagePixelType tetr_label,ImagePixelType linked_tetr_label)
 {
	// enum TRI_TYPE {EXTERIOR_SURFACE = 0, INTERIOR_SURFACE = 1, INTERIOR = 2};
	if(triType != 0 && triType != 1)
		itkExceptionMacro("Invalid type of triangle found!");

	TetrahedronType* pTestTetr = NULL;

	// EXTERIOR_SURFACE = 0
	if(triType == 0)
	{
		// The orientation of the surface triangle is always checked
		// with the unique tet that the triangle belongs to.
		if(!pTetr)
			itkExceptionMacro("Tetr is NULL!");

		pTestTetr = pTetr;
	}
	// INTERIOR_SURFACE = 1
	else
	{
		// The orientation of the surface triangle is always checked
		// with the tet with the smaller label that the triangle belongs to
		// (the interior surface triangle is always common between two tets).
		if(!pTetr || !pLinkedTetr)
			itkExceptionMacro("Tetr or pLinkedTetr is NULL!");

		if(tetr_label == linked_tetr_label)
			itkExceptionMacro("The tatrahedral labels are the same in an interior surface triangle!");

		// Choose the smaller label
		if(tetr_label < linked_tetr_label)
			pTestTetr = pTetr;
		else
			pTestTetr = pLinkedTetr;
	}

	unsigned int fourthVertexId=0;
	if(GetFourthVertexId(pTestTetr,ids,fourthVertexId) == false)
		itkExceptionMacro("Triangle's ids do not exist in the Tetrahedron!");

	PointType pos0 = this->GetInput()->GetPoint(ids[0]);
	PointType pos1 = this->GetInput()->GetPoint(ids[1]);
	PointType pos2 = this->GetInput()->GetPoint(ids[2]);
	PointType pos3 = this->GetInput()->GetPoint(fourthVertexId);

	VectorType vec0,vec1,vec2,vec3;
	for(unsigned int i=0; i<MeshDimension; i++)
	{
		vec0[i] = pos0[i];
		vec1[i] = pos1[i];
		vec2[i] = pos2[i];
		vec3[i] = pos3[i];
	}

	double volume = CalculateTetrahedronVolume(vec0.GetVnlVector(),vec1.GetVnlVector(),vec2.GetVnlVector(),vec3.GetVnlVector());
	if(volume < 0.0)
	{
		// Rearrange the order of the vertices ids so that : volume > 0
		unsigned int tempVal = ids[1];
		ids[1] = ids[2];
		ids[2] = tempVal;
	}

	// check again the volume for verification
	pos0 = this->GetInput()->GetPoint(ids[0]);
	pos1 = this->GetInput()->GetPoint(ids[1]);
	pos2 = this->GetInput()->GetPoint(ids[2]);
	pos3 = this->GetInput()->GetPoint(fourthVertexId);
	for(unsigned int i=0; i<MeshDimension; i++)
	{
		vec0[i] = pos0[i];
		vec1[i] = pos1[i];
		vec2[i] = pos2[i];
		vec3[i] = pos3[i];
	}

	volume = CalculateTetrahedronVolume(vec0.GetVnlVector(),vec1.GetVnlVector(),vec2.GetVnlVector(),vec3.GetVnlVector());
	if(volume < 0.0)
		itkExceptionMacro("Check Surface Triangle Orientation failed!");
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::BuildSourcePointsConnectivity()
 {
	std::cout << "Building Surface Vertices Connectivity..." << std::endl;

	this->m_SurfaceVertexToNeighborVerticesMapContainer = SurfaceVertexToNeighborVerticesMapContainerType::New();
	SurfaceVertexToNeighborVerticesMapContainerType::iterator mapIter;
	typename SurfaceTrianglesContainerType::Iterator triIter;

	std::set<unsigned int>* pNeighborVerticesSet;
	VectorTriIdType vTriId;

	for(triIter = this->m_SurfaceTrianglesContainer->Begin(); triIter != this->m_SurfaceTrianglesContainer->End(); triIter++)
	{
		vTriId = triIter.Value();

		unsigned int id;
		unsigned int idOther1,idOther2;

		for(unsigned int i=0; i<MeshDimension; i++)
		{
			id = vTriId[i];
			mapIter = this->m_SurfaceVertexToNeighborVerticesMapContainer->CastToSTLContainer().find(id);

			// Update the source points connectivity
			if(mapIter != this->m_SurfaceVertexToNeighborVerticesMapContainer->CastToSTLContainer().end())
			{
				pNeighborVerticesSet = mapIter->second;
			}
			// Create new connectivity
			else
			{
				pNeighborVerticesSet = new std::set<unsigned int>();
				this->m_SurfaceVertexToNeighborVerticesMapContainer->CastToSTLContainer()[id] = pNeighborVerticesSet;
			}

			GetOtherVerticesId(vTriId,id,idOther1,idOther2);

			// Insert the other two vertices of the triangle
			pNeighborVerticesSet->insert(idOther1);
			pNeighborVerticesSet->insert(idOther2);
		}
	}
	std::cout << "Size 0f m_Vertex2NeighborVerticesMap : " << m_SurfaceVertexToNeighborVerticesMapContainer->Size() << std::endl;
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::PrintTargetPointsITK(typename TargetPointSetType::Pointer pPointSet,ImageType* pImage, std::string fname)
 {
	std::cout << "Writing Target Points Image at : "<< fname << std::endl;
	if(pPointSet.IsNull())
		return;

	typename TargetPointSetType::PointsContainer* pPointsContainer = pPointSet->GetPoints();
	if(!pPointsContainer || pPointsContainer->Size() == 0 || !pImage || fname.empty())
		return;

	typename ImageType::Pointer pOutputImage = ImageType::New();

	pOutputImage->SetBufferedRegion( pImage->GetBufferedRegion() );
	pOutputImage->SetRequestedRegion( pImage->GetRequestedRegion() );
	pOutputImage->SetLargestPossibleRegion( pImage->GetLargestPossibleRegion() );
	pOutputImage->SetDirection(pImage->GetDirection());
	pOutputImage->SetOrigin(pImage->GetOrigin());
	pOutputImage->SetSpacing(pImage->GetSpacing());
	pOutputImage->Allocate();
	pOutputImage->FillBuffer(0);
	pOutputImage->Update();

	typedef TargetPointSetIteratorType PointIteratorType;
	PointIteratorType pointItr = pPointsContainer->Begin();
	PointIteratorType pointItrEnd = pPointsContainer->End();
	ImageIndexType index;

	unsigned int NumVoxels = 0;
	while ( pointItr != pointItrEnd )
	{
		if ( pOutputImage->TransformPhysicalPointToIndex(pointItr.Value(), index) )
		{
			pOutputImage->SetPixel(index,1);
			NumVoxels++;
		}
		else
		{
			std::cout << "point : " << pointItr.Value() << std::endl;
			std::cout << "index : " << index << std::endl;
			itkExceptionMacro("A point is outside the image!");
		}
		pointItr++;
	}

	std::cout << "Num Voxels: "<< NumVoxels << std::endl;
	WriteImage<ImagePixelType,MeshDimension>(pOutputImage,fname);
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::PrintTargetPointsVTK(typename TargetPointSetType::Pointer pTargetPointSet,std::string fName)
 {
	std::cout << "Printing Target Points..." << std::endl;
	if(!pTargetPointSet || fName.empty())
	{
		std::cout << "Invalid Input!" << std::endl;
		return;
	}

	vtkNew<vtkUnstructuredGrid> pUG;
	vtkPoints* pPoints = vtkPoints::New();

	// Set coordinates
	typename TargetPointSetType::PointsContainer* pPointsContainer = pTargetPointSet->GetPoints();
	if(!pPointsContainer)
	{
		std::cout << "Points Container is NULL!" << std::endl;
		return;
	}

	int numPoints = pPointsContainer->Size();
	typename TargetPointSetType::PointType pt;

    for(unsigned int i=0; i < (unsigned int)numPoints; i++ )
    {
    	pt = pPointsContainer->ElementAt(i);
    	pPoints->InsertPoint(i,pt[0],pt[1],pt[2]);
    }

    // Set the points to Unstructured Grid
	pUG->SetPoints(pPoints);
	pPoints->Delete();

	vtkIdList* ptIds = vtkIdList::New();

    for(unsigned int i = 0; i < (unsigned int)numPoints; i++)
	{
		ptIds->InsertId(0,i);

		// Insert the cell to Unstructured Grid
		pUG->InsertNextCell(1,ptIds);
	}
	ptIds->Delete();

	vtkUnstructuredGridWriter* pWriter = vtkUnstructuredGridWriter::New();
	pWriter->SetFileTypeToASCII();
	pWriter->SetFileName(fName.c_str());
	pWriter->SetInputData(pUG.GetPointer());
	pWriter->Write();
	pWriter->Update();

	pWriter->Delete();
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::PrintSourcePointsVTK(typename SurfaceVertexToLabelMapContainerType::Pointer pSurfVerticesMapContainer,const typename MeshType::PointsContainer* pPointsContainer,std::string fName)
 {
	std::cout << "Printing Source Points..." << std::endl;

	if(!pSurfVerticesMapContainer)
		return;
	if(!pPointsContainer)
		return;

	vtkNew<vtkUnstructuredGrid> pUG;
	vtkPoints* pPoints = vtkPoints::New();
	int numPoints = pPointsContainer->Size();
	PointType pt;
	for( int i=0; i < numPoints; i++ )
	{
		pt = pPointsContainer->ElementAt(i);
		pPoints->InsertPoint(i,pt[0],pt[1],pt[2]);
	}
    pUG->SetPoints(pPoints);
    pPoints->Delete();

	vtkIdList* ptIds = vtkIdList::New();
	long int id;

	typename SurfaceVertexToLabelMapContainerType::iterator it;
	for(it = pSurfVerticesMapContainer->CastToSTLContainer().begin(); it != pSurfVerticesMapContainer->CastToSTLContainer().end(); it++)
	{
		id = it->first;
		ptIds->InsertId(0,id);

		// Insert the cell to Unstructured Grid
		pUG->InsertNextCell(1,ptIds);
	}
	ptIds->Delete();


	vtkUnstructuredGridWriter* pWriter = vtkUnstructuredGridWriter::New();
	pWriter->SetFileTypeToASCII();
	pWriter->SetFileName(fName.c_str());
	pWriter->SetInputData(pUG.GetPointer());
	pWriter->Write();
	pWriter->Update();
	pWriter->Delete();
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::PrintSurfaceTriangles(typename SurfaceTrianglesContainerType::Pointer pTrianglesContainer,const typename MeshType::PointsContainer* pPointsContainer,std::string fName)
 {
	std::cout << "Printing Surface Triangles..." << std::endl;
	if(!pTrianglesContainer)
		return;
	if(!pPointsContainer)
		return;

	vtkNew<vtkUnstructuredGrid> pUG;
	vtkPoints* pPoints = vtkPoints::New();
	int numPoints = pPointsContainer->Size();
	PointType pt;
	for( int i=0; i < numPoints; i++ )
	{
		pt = pPointsContainer->ElementAt(i);
		pPoints->InsertPoint(i,pt[0],pt[1],pt[2]);
	}
    pUG->SetPoints(pPoints);
    pPoints->Delete();

	vtkIdList* ptIds = vtkIdList::New();

	VectorTriIdType vTriId;
	typename SurfaceTrianglesContainerType::Iterator it;
	for(it = pTrianglesContainer->Begin(); it != pTrianglesContainer->End(); it++)
	{
		vTriId = it.Value();

		ptIds->InsertId(0,vTriId[0]);
		ptIds->InsertId(1,vTriId[1]);
		ptIds->InsertId(2,vTriId[2]);

		// Insert the cell to Unstructured Grid
		pUG->InsertNextCell(5,ptIds);
	}
	ptIds->Delete();

	vtkUnstructuredGridWriter* pWriter = vtkUnstructuredGridWriter::New();
	pWriter->SetFileTypeToASCII();
	pWriter->SetFileName(fName.c_str());
	pWriter->SetInputData(pUG.GetPointer());
	pWriter->Write();
	pWriter->Update();
	pWriter->Delete();
 }


template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::InitializeMaterials()
 {
	MaterialContainerType * materialContainer = this->GetFEMObject()->GetMaterialContainer();

	if(!materialContainer)
		itkExceptionMacro("Missing material container");

	materialContainer->Initialize();

	// fix material to linear elasticity
	this->SetMaterials();
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::SetMaterials()
 {
	if(!this->GetFEMObject())
		return;

	// In the future we will support N materials.
	this->m_MaterialBrain = MaterialType::New();

	// Properties taken from Miga's paper : Modeling of retraction and resection fot intraoperative Updating images
	if(this->GetYoungModulus() < 1e-06)
		itkExceptionMacro("Invalid value for Young's Modulus : " << this->GetYoungModulus());

	if(this->GetPoissonRatio() < 1e-06 || this->GetPoissonRatio() > 0.499)
		itkExceptionMacro("Invalid value for Poisson's Ratio: " << this->GetPoissonRatio());

	this->GetMaterialBrain()->SetYoungsModulus(this->GetYoungModulus());
	this->GetMaterialBrain()->SetPoissonsRatio(this->GetPoissonRatio());
	this->GetFEMObject()->AddNextMaterial(this->m_MaterialBrain);
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::InitializeFEMObject()
 {
	std::cout << "\nInitializing FEMObject..." << std::endl;

	if(this->GetFEMObject())
		return;

	this->m_FEMObject = FEMObjectType::New();
	this->m_FEMObject->Initialize();

	this->InitializeMaterials();
	this->InitializeNodes();
	this->InitializeElements();
	this->InitializeLoads();

	// produce DOF
	this->GetFEMObject()->FinalizeMesh();
 }

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::InitializeNodes()
 {
	NodeContainerType *nodeContainer = this->GetFEMObject()->GetNodeContainer();

	if(!nodeContainer)
		itkExceptionMacro("Missing node container");

	nodeContainer->Initialize();
	const PointsContainerType* meshPoints = this->GetInput()->GetPoints();
	PointsIteratorType it = meshPoints->Begin();
	FEMVectorType point(MeshDimension);

	// Initialize nodes
	while(it != meshPoints->End())
	{
		for(unsigned i = 0; i < MeshDimension; i++)
			point[i]= it.Value()[i];

		NodeType::Pointer node = NodeType::New();
		node->SetCoordinates(point);
		node->SetGlobalNumber(it.Index());

		this->GetFEMObject()->AddNextNode(node);
		++it;
	}
	if( this->GetFEMObject()->GetNumberOfNodes() != this->GetInput()->GetPoints()->Size())
		itkExceptionMacro("FEMObject generated invalid number of Nodes! ");

	std::cout << "Num Nodes at FEMObject : " << this->GetFEMObject()->GetNumberOfNodes() << std::endl;
 }


template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::InitializeElements()
 {
	ElementContainerType *elementContainer = this->GetFEMObject()->GetElementContainer();
	if(!elementContainer)
		itkExceptionMacro("Missing element container");
	elementContainer->Initialize();

	const MeshType* pMesh = this->GetInput();
	if(!pMesh)
		itkExceptionMacro("The input mesh is NULL!");

	this->m_TetrahedraDataContainer = TetrahedraDataContainerType::New();
	this->m_TetrahedraDataContainer->Initialize();

	//Mesh cell iterator
	CellIteratorType cellIterator = pMesh->GetCells()->Begin();
	CellIteratorType cellEnd = pMesh->GetCells()->End();

	if(this->GetInput()->GetCellData() == NULL)
		itkExceptionMacro("Undefined cell data for the itk::Mesh");

	if(this->GetInput()->GetCellData()->Size() != pMesh->GetCells()->Size())
		itkExceptionMacro("Different size between cell container and cells data container for the itk::Mesh");

	CellDataIteratorType CellDataIterator = this->GetInput()->GetCellData()->Begin();
	CellType* cell = NULL;
	typename MeshType::CellIdentifier cellId = 0;
	ImagePixelType label;

	// Loop over all tets
	while( cellIterator != cellEnd )
	{
		cell = cellIterator.Value();
		label = CellDataIterator.Value();

		// Only tetrahedra
		if(cell->GetType() != CellType::TETRAHEDRON_CELL)
			itkExceptionMacro("Only Tetrahedral cells are supported. The current element type is : " << cell->GetType());

		FEMTetrahedronType::Pointer tetrahedronEle = FEMTetrahedronType::New();
		TetrahedronType * pTetr = static_cast<TetrahedronType *>( cell );

		PointIdIterator pointIdIter = pTetr->PointIdsBegin();
		PointIdIterator pointIdEnd = pTetr->PointIdsEnd();
		NodeType::Pointer pNode;
		NodeType::SetOfElements adjacentElements;
		typename MeshType::CellFeatureIdentifier pointId = 0;
		while( pointIdIter != pointIdEnd )
		{
			pNode = this->GetFEMObject()->GetNode(*pointIdIter);
			// set the node for the tetrahedron
			tetrahedronEle->SetNode(pointId,pNode);

			// Update the list of pointers to elements that use this node
			adjacentElements = pNode->m_elements;
			adjacentElements.insert(tetrahedronEle);
			pNode->m_elements = adjacentElements;

			++pointIdIter;
			pointId++;
		}

		((Element::Pointer)tetrahedronEle)->PopulateEdgeIds();
		tetrahedronEle->SetGlobalNumber(cellId);
		tetrahedronEle->SetMaterial( static_cast<MaterialType *>( this->GetFEMObject()->GetMaterial(0).GetPointer() ) );
		this->GetFEMObject()->AddNextElement(tetrahedronEle.GetPointer());
		this->m_TetrahedraDataContainer->InsertElement(cellId,label);

		cellId++;
		++cellIterator;
		++CellDataIterator;
	}
	if( this->GetFEMObject()->GetNumberOfElements() != pMesh->GetCells()->Size())
		itkExceptionMacro("FEMObject generated invalid number of Tets! ");

	std::cout << "Num Tets at FEMObject : " << this->GetFEMObject()->GetNumberOfElements() << std::endl;
}

template<class TInputMesh, class TInputImage, class TOutputMesh>
void FEMCompressor<TInputMesh, TInputImage, TOutputMesh>
::InitializeLoads()
 {
	LoadContainerType *loadContainer = this->GetFEMObject()->GetLoadContainer();
	if(!loadContainer)
		itkExceptionMacro("Missing load container");
	loadContainer->Initialize();

	unsigned int id = 0;
	Element::VectorType shape;

	NodeType::Pointer pNode;
	NodeType::SetOfElements::iterator it;
	Element* pElem = NULL;
	typename SurfaceVertexToLabelMapContainerType::Iterator iter;
	unsigned int globalNum = 0;

	// Loop over all surface mesh vertices
	for(iter = this->m_SurfaceVertexToLabelMapContainer->Begin(); iter != this->m_SurfaceVertexToLabelMapContainer->End(); iter++)
	{
		id = iter.Index();
		pNode = this->GetFEMObject()->GetNode(id);

		LoadType::Pointer pLoad =  LoadType::New();
		pLoad->SetGlobalNumber(globalNum);
		pLoad->SetConfidence(1); // confidence 1
		pLoad->SetIsOutOfMesh(false);
		pLoad->SetOutlier(false);

		// As contained element we set the first of the adjucent elements of the vertex
		it = pNode->m_elements.begin();
		pElem = *it;
		if(!pElem)
			itkExceptionMacro("The adjacent element of a node is NULL!");

		// We set exact (0 or 1) the shape functions because we know that the point belongs to a vertex!
		// In that way we avoid numerical errors
		shape.set_size(pElem->GetNumberOfNodes());
		unsigned int hits=0;
		for(unsigned int i=0; i<pElem->GetNumberOfNodes(); i++)
		{
			if(id != (unsigned int)pElem->GetNode(i)->GetGlobalNumber())
				shape[i] = 0.0;
			else
			{
				shape[i] = 1.0;
				hits++;
			}
		}
		// Check the shape
		if(hits != 1)
			itkExceptionMacro("Invalide shape function!");
		//std::cout << shape << std::endl;

		pLoad->SetShape(shape);
		pLoad->SetContainedElement(pElem);
		this->GetFEMObject()->AddNextLoad(fem::Load::Pointer(pLoad));
		globalNum++;
	}
	std::cout << "Num Loads at FEMObject : " << this->GetFEMObject()->GetNumberOfLoads() << std::endl;
 }








} // fem
} // itk


#endif /* ITKFEMCOMPRESSOR_HXX_ */
