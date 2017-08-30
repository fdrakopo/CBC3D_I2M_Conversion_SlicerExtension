/*
 * itkFEMCompressor.h
 *
 */

#ifndef ITKFEMCOMPRESSOR_H_
#define ITKFEMCOMPRESSOR_H_

#include "vnl/vnl_cross.h"

#include "itkMeshToMeshFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkFEMObject.h"
#include "itkFEMElement3DC0LinearTetrahedronStrain.h"
#include "itkTetrahedronCell.h"
#include "itkFEMCompressorSolver.h"
#include "itkFEMLoadNoisyLandmark.h"

namespace itk
{ 
namespace fem
{

template<class TInputMesh, class TInputImage, class TOutputMesh>
class ITK_EXPORT FEMCompressor : public MeshToMeshFilter<TInputMesh,TOutputMesh>
{

public:
	typedef FEMCompressor                						  	Self;
	typedef MeshToMeshFilter<TInputMesh,TOutputMesh>   				Superclass;
	typedef SmartPointer<Self>                                    	Pointer;
	typedef SmartPointer<const Self>                              	ConstPointer;


	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Extract dimension from the output mesh. */
	itkStaticConstMacro( MeshDimension, unsigned int, TOutputMesh::PointDimension);

	typedef TInputMesh                    							MeshType;
	typedef typename MeshType::PixelType       						MeshPixelType;
	typedef typename MeshType::Pointer         						MeshPointerType;
	typedef typename MeshType::CellType                           	CellType;
	typedef typename CellType::CellAutoPointer                   	CellAutoPointerType;
	typedef typename MeshType::CellsContainer                    	CellsContainerType;
	typedef typename CellsContainerType::ConstIterator              CellIteratorType;
	typedef typename MeshType::CellDataContainer::ConstIterator 	CellDataIteratorType;

	typedef TInputImage               								ImageType;
	typedef typename ImageType::Pointer         					ImagePointerType;
	typedef typename ImageType::PixelType         					ImagePixelType;
	typedef typename ImageType::IndexType         					ImageIndexType;
	typedef typename ImageType::SpacingType         				ImageSpacingType;
	typedef typename ImageType::RegionType 	        				ImageRegionType;
	typedef typename ImageType::SizeType 	        				ImageSizeType;
	typedef typename std::set<ImagePixelType>						PixelTypeSet;

	typedef TOutputMesh                    							OutputMeshType;
	typedef typename OutputMeshType::Pointer         				OutputMeshPointerType;

	typedef typename itk::PointSet<PixelTypeSet>					TargetPointSetType;
	typedef typename TargetPointSetType::PointsContainerIterator	TargetPointSetIteratorType;
	typedef typename TargetPointSetType::PointDataContainerIterator	TargetPointDataSetIteratorType;

	typedef typename MeshType::PointType 							PointType;
	typedef typename MeshType::PointsContainer						PointsContainerType;
	typedef typename PointsContainerType::ConstIterator				PointsIteratorType;
//	typedef typename PointSetType::PointDataContainer				PointsDataContainerType;
//	typedef typename PointsDataContainerType::ConstIterator			PointDataIteratorType;

	typedef Image< unsigned char, MeshDimension>  					SelectionMapImageType;
	typedef Offset< MeshDimension >       							OffsetType;


	typedef FEMObject<MeshDimension>                 				FEMObjectType;
	typedef typename FEMObjectType::Pointer							FEMObjectPointerType;
	typedef LoadNoisyLandmark 										LoadType;
	typedef CompressorSolver<ImagePixelType,MeshDimension>			FEMCompressorSolverType;
	typedef typename FEMCompressorSolverType::Pointer				FEMCompressorSolverPointerType;


	/** FEM container typedef support */
	typedef typename FEMObjectType::LoadContainerType        		LoadContainerType;
	typedef typename FEMObjectType::NodeContainerType        		NodeContainerType;
	typedef typename FEMObjectType::ElementContainerType     		ElementContainerType;
	typedef typename FEMObjectType::MaterialContainerType    		MaterialContainerType;

	/** FEM element typedef support */
	typedef Element::VectorType                       				FEMVectorType;

	/** FEM element typedef support */
	typedef Element3DC0LinearTetrahedronStrain        				FEMTetrahedronType;
	typedef TetrahedronCell<CellType>                         		TetrahedronType;
	typedef TriangleCell<CellType>		                       		TriangleType;

	typedef typename CellType::PointIdIterator                		PointIdIterator;

	/** FEM node typedef support */
	typedef Element::Node                             				NodeType;

	/** FEM material typedef support */
	typedef MaterialLinearElasticity                  				MaterialType;
	typedef MaterialType::Pointer                     				MaterialPointerType;


	typedef typename itk::Vector<double,MeshDimension>												VectorType;
	typedef typename itk::Vector<unsigned int,MeshDimension>										VectorTriIdType;
	typedef typename itk::MapContainer<unsigned int, PixelTypeSet>								    SurfaceVertexToLabelMapContainerType;
	typedef typename itk::MapContainer<unsigned int, std::set<unsigned int>* > 						SurfaceVertexToNeighborVerticesMapContainerType;
	typedef typename itk::VectorContainer<unsigned int,	VectorTriIdType >							SurfaceTrianglesContainerType;
	typedef typename itk::VectorContainer<unsigned int,	ImagePixelType>								TetrahedraDataContainerType;

	itkSetObjectMacro(Image, ImageType);
	itkGetObjectMacro(Image, ImageType);

	itkSetMacro(Lambda, float);
	itkGetMacro(Lambda, float);

	itkSetMacro(YoungModulus, float);
	itkGetMacro(YoungModulus, float);

	itkSetMacro(PoissonRatio, float);
	itkGetMacro(PoissonRatio, float);

	itkSetMacro(NumIterations, unsigned int);
	itkGetMacro(NumIterations, unsigned int);

	itkSetObjectMacro(FEMObject, FEMObjectType);
	itkGetObjectMacro(FEMObject, FEMObjectType);

	itkSetObjectMacro(MaterialBrain, MaterialType);
	itkGetObjectMacro(MaterialBrain, MaterialType);

	// set/get non-connectivity with radius = 1 of dimension connect.
	// 0 is vertex non-connectivity (26 in 3D)
	// 1 is edge non-connectivity (18 in 3D)
	// 2 is face non-connectivity (6 in 3D)
	// 3 is NO non-connectivity */
	itkSetMacro(NonConnectivityType, unsigned int);
	itkGetMacro(NonConnectivityType, unsigned int);

	itkSetMacro(OutputSourcePointsFileName, std::string);
	itkGetMacro(OutputSourcePointsFileName, std::string);

	itkSetMacro(OutputTargetPointsFileName, std::string);
	itkGetMacro(OutputTargetPointsFileName, std::string);

	itkSetMacro(OutputSurfaceTrianglesFileName, std::string);
	itkGetMacro(OutputSurfaceTrianglesFileName, std::string);

	void GenerateData();

protected:

	FEMCompressor();
	~FEMCompressor();

	/** connectivity constants */
	enum
	{
		VERTEX_NON_CONNECTIVITY = 0,
		EDGE_NON_CONNECTIVITY = 1,
		FACE_NON_CONNECTIVITY = 2,
		NO_NON_CONNECTIVITY = 3
	};

	void ExtractSourcePoints();
	void BuildSourcePointsConnectivity();
	void ExtractTargetPoints();
	void CheckSurfaceTriangleOrientation(int triType,unsigned int ids[3],TetrahedronType* pTetr,TetrahedronType* pLinkedTetr,ImagePixelType tetr_label,ImagePixelType linked_tetr_label);

	void InitializeFEMObject();
	void InitializeMaterials();
	void InitializeNodes();
	void InitializeElements();
	void InitializeLoads();
	void SetMaterials();
	void PrintSurfaceTriangles(typename SurfaceTrianglesContainerType::Pointer pTrianglesContainer,const typename MeshType::PointsContainer* pPointsContainer,std::string fName);
	void PrintSourcePointsVTK(typename SurfaceVertexToLabelMapContainerType::Pointer pSurfVerticesMapContainer,const typename MeshType::PointsContainer* pPointsContainer,std::string fName);
	void PrintTargetPointsVTK(typename TargetPointSetType::Pointer pTargetPointSet,std::string fName);
	void PrintTargetPointsITK(typename TargetPointSetType::Pointer pPointSet,ImageType* pImage, std::string fname);
	void UpdateOutputMesh();
	void MarkOffConnectedPoints(ImageIndexType index);
	void ComputeConnectivityOffsets();

	bool SameTriangles(VectorTriIdType vTri1Id,VectorTriIdType vTri2Id)
	{
		unsigned int hits = 0;
		for(unsigned int i=0; i<MeshDimension; i++)
		{
			for(unsigned int j=0; j<MeshDimension; j++)
			{
				if(vTri1Id[i] == vTri2Id[j])
					hits++;
			}
		}
		if(hits > 3)
			itkExceptionMacro("Invaled input triangles ids!");

		if(hits == 3)
			return true;
		else
			return false;
	}

	bool SameTriangles(TriangleType* pTriA,TriangleType* pTriB)
	{
		if(!pTriA || !pTriB)
			return false;

		if(pTriA == pTriB)
			return true;

		PointIdIterator IterA = pTriA->PointIdsBegin();
		PointIdIterator EndA =  pTriA->PointIdsEnd();
		PointIdIterator IterB = pTriB->PointIdsBegin();
		PointIdIterator EndB =  pTriB->PointIdsEnd();

		unsigned int hits = 0;
		while( IterA != EndA )
		{
			while( IterB != EndB )
			{
				if(*IterA == *IterB)
				{
					hits++;
					break;
				}

				IterB++;
			}
			IterA++;
			IterB = pTriB->PointIdsBegin();
		}
		if(hits == 3)
			return true;
		else
			return false;
	}
	void CheckForDuplicateSurfaceTriangles()
	{
		std::cout << "Checking For Duplicate Surface Triangles..." << std::endl;
		if(this->m_SurfaceTrianglesContainer)
		{
			for(int i=0; i< this->m_SurfaceTrianglesContainer->Size(); i++)
			{
				VectorTriIdType vTri1Id = this->m_SurfaceTrianglesContainer->GetElement(i);
				for(int j=i+1; j< this->m_SurfaceTrianglesContainer->Size(); j++)
				{
					VectorTriIdType vTri2Id = this->m_SurfaceTrianglesContainer->GetElement(j);
					if(SameTriangles(vTri1Id,vTri2Id))
						itkExceptionMacro("Same triangles were found!");
				}
			}
		}
	}
	void CheckForDuplicateSurfaceVertices()
	{
		std::cout << "Checking For Duplicate Surface Vertices..." << std::endl;
		if(this->m_SurfaceVertexToLabelMapContainer)
		{
			typename SurfaceVertexToLabelMapContainerType::iterator it1;
			typename SurfaceVertexToLabelMapContainerType::iterator it2;
			int count=1;
			for(it1 = this->m_SurfaceVertexToLabelMapContainer->CastToSTLContainer().begin(); it1 != this->m_SurfaceVertexToLabelMapContainer->CastToSTLContainer().end(); it1++)
			{
				unsigned int id1 = it1->first;

				it2 = this->m_SurfaceVertexToLabelMapContainer->CastToSTLContainer().begin();
				std::advance (it2,count);

				while(it2 != this->m_SurfaceVertexToLabelMapContainer->CastToSTLContainer().end())
				{
					unsigned int id2 = it2->first;
					if(id1 == id2)
						itkExceptionMacro("Same Surface vertices were found!");
					it2++;
				}
				count++;
			}
		}
	}
	void GetOtherVerticesId(VectorTriIdType vTriId,unsigned int id,unsigned int& idOther1,unsigned int& idOther2)
	{
		if(id == vTriId[0])
		{
			idOther1 = vTriId[1];
			idOther2 = vTriId[2];
		}
		else if(id == vTriId[1])
		{
			idOther1 = vTriId[2];
			idOther2 = vTriId[0];
		}
		else if(id == vTriId[2])
		{
			idOther1 = vTriId[0];
			idOther2 = vTriId[1];
		}
		else
		{
			itkExceptionMacro("Invalid vertex id's!");
		}
	}

	bool GetFourthVertexId(TetrahedronType* pTetr,unsigned int ids[3],unsigned int& idOther)
	{
		if(pTetr->GetNumberOfPoints() != 4)
			itkExceptionMacro("Invalid number of points for tetrahedron cell!");

		PointIdIterator pIter = pTetr->PointIdsBegin();
		PointIdIterator pIterEnd = pTetr->PointIdsEnd();

		unsigned int hits = 0;

		while(pIter != pIterEnd)
		{
			if( *pIter == ids[0] || *pIter == ids[1] || *pIter == ids[2] )
				hits++;
			else
				idOther = *pIter;

			pIter++;
		}

		if(hits != 3)
			return false;
		else
			return true;
	}

	bool SameLabelSets(PixelTypeSet set1,PixelTypeSet set2)
	{
		if(set1.size() != set2.size())
			return false;

		typename PixelTypeSet::iterator it1 = set1.begin();
		typename PixelTypeSet::iterator it2 = set2.begin();

		while(it1 != set1.end())
		{
			if(*it1 != *it2)
				return false;

			++it1;
			++it2;
		}
		return true;
	}

	unsigned int GetNumOfCommonIds(unsigned int vTriId[3],unsigned int vTetId[4])
	{
		unsigned int num = 0;
		for(unsigned int i=0; i<3; i++)
		{
			if(vTriId[i] == vTetId[0] ||
			   vTriId[i] == vTetId[1] ||
			   vTriId[i] == vTetId[2] ||
			   vTriId[i] == vTetId[3])
				num++;
		}
		return num;
	}

	template<class TPixelType>
	bool CheckTetrahedraVolumes(const typename itk::Mesh<TPixelType,MeshDimension>* pMesh)
	{
		const double ZERO6 = 0.000001;

		std::cout << "Checking tetrahedra Volumes..." << std::endl;
		if(!pMesh)
			return false;

		typedef itk::Mesh<TPixelType,MeshDimension>				MeshType;
		typedef typename MeshType::CellType             CellType;
		typedef typename CellType::PointIdIterator    	PointIdIterator;


		//Mesh cell iterator
		typedef typename MeshType::CellsContainer::ConstIterator     CellIterator;
		CellIterator cellIterator = pMesh->GetCells()->Begin();
		CellIterator cellEnd = pMesh->GetCells()->End();
		vnl_vector<double> pos0(3),pos1(3),pos2(3),pos3(3);
		vnl_vector<double> pos(3);
		typename itk::fem::FEMObject<MeshDimension>::NodeIdentifier nodeId;
		itk::Point<double,MeshDimension> point;

		double vol = 0;
		unsigned int invertedTetrahedra = 0;
		unsigned int count=0;

		while( cellIterator != cellEnd )
		{
			CellType * cell = cellIterator.Value();
			if(cell->GetType() != CellType::TETRAHEDRON_CELL)
			{
				std::cout << "ERROR at CheckTetrahedraVolumes : Cell type is not tetrahedron!" << std::endl;
				return false;
			}

			PointIdIterator pointIdIter = cell->PointIdsBegin();
			PointIdIterator pointIdEnd = cell->PointIdsEnd();

			count=0;
			while( pointIdIter != pointIdEnd )
			{
				nodeId = *pointIdIter;
				point = pMesh->GetPoint(nodeId);

				for(unsigned int i=0; i<MeshDimension; i++)
					pos[i] = point[i];

				if(count==0)
					pos0 = pos;
				else if(count==1)
					pos1 = pos;
				else if(count==2)
					pos2 = pos;
				else
					pos3 = pos;

				count++;
				++pointIdIter;
			}

			vol = CalculateTetrahedronVolume(pos0,pos1,pos2,pos3);
			if(vol < ZERO6)
				invertedTetrahedra++;

			//cout << "Volume " << vol << endl;

			++cellIterator;
		}

		if(invertedTetrahedra == 0)
		{
			//		std::cout << "OK" << std::endl;
			return true;
		}
		else
		{
			std::cout <<  invertedTetrahedra << " inverted tetrahedra were found!" << std::endl;
			return false;
		}
	}

	double CalculateTetrahedronVolume(vnl_vector<double> pos0,vnl_vector<double> pos1,vnl_vector<double> pos2,vnl_vector<double> pos3)
	{
		vnl_vector<double> p = pos1 - pos0;
		vnl_vector<double> q = pos2 - pos0;
		vnl_vector<double> r = pos3 - pos0;

		double volume = ( dot_product(vnl_cross_3d(p,q),r) ) /6.0;
		return volume;
	}

	template<class TPixel, unsigned int  VImageDimension >
	void WriteImage(typename itk::Image< TPixel, VImageDimension >::Pointer pImage ,std::string fname)
	{
		printf("\nSaving Image at : %s\n",fname.c_str());
		typedef itk::ImageFileWriter< itk::Image< TPixel,VImageDimension > > WriterType;
		typename WriterType::Pointer pWriter = WriterType::New();
		pWriter->SetFileName(fname);
		pWriter->SetInput(pImage);
		pWriter->Update();
	}


private:

	ImagePointerType				m_Image;
	unsigned int          			m_NonConnectivityType;
    std::vector< OffsetType >  		m_NonConnectivityOffsets;
	float 							m_Lambda;
	unsigned int 					m_NumIterations;
	FEMObjectPointerType			m_FEMObject;
	MaterialPointerType 			m_MaterialBrain;
	float 							m_YoungModulus;
	float 							m_PoissonRatio;
	FEMCompressorSolverPointerType	m_FEMCompressorSolver;

	typename SurfaceVertexToLabelMapContainerType::Pointer				m_SurfaceVertexToLabelMapContainer;
	typename SurfaceVertexToNeighborVerticesMapContainerType::Pointer 	m_SurfaceVertexToNeighborVerticesMapContainer;
	typename SurfaceTrianglesContainerType::Pointer						m_SurfaceTrianglesContainer;
	typename TargetPointSetType::Pointer								m_TargetPointSet;
	typename SelectionMapImageType::Pointer								m_SelectionMapImage;

	// Helpers for printing
	typename TetrahedraDataContainerType::Pointer						m_TetrahedraDataContainer;
	std::string 														m_OutputSourcePointsFileName;
	std::string 														m_OutputTargetPointsFileName;
	std::string 														m_OutputSurfaceTrianglesFileName;
};


} // fem
} // itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFEMCompressor.hxx"
#endif


#endif /* ITKFEMCOMPRESSOR_H_ */
