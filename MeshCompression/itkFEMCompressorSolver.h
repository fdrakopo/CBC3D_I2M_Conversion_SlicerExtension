/*
 * itkFEMCompressorSolver.h
 *
 */

#ifndef ITKFEMCOMPRESSORSOLVER_H_
#define ITKFEMCOMPRESSORSOLVER_H_


#include "itkFEMSolver.h"
#include "itkFEMElementBase.h"
#include "itkFEMMaterialBase.h"
#include "itkFEMLoadBase.h"
#include "itkFEMLinearSystemWrapperVNL.h"
#include "itkFEMLinearSystemWrapperItpack.h"
#include "itkFEMLoadNoisyLandmark.h"


namespace itk
{
namespace fem
{

template <class TPixel,unsigned int VDimension>
class ITK_EXPORT CompressorSolver : public Solver<VDimension>
{
public:
  /** Standard class typedefs. */
  typedef CompressorSolver          Self;
  typedef Solver<VDimension>        Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(CompressorSolver, Solver);

  itkStaticConstMacro(FEMDimension, unsigned int, VDimension);

  /** Inherit some types from the superclass. */
  typedef typename Superclass::VectorType                      	VectorType;
  typedef typename Superclass::Float                           	Float;

  typedef TPixel												ImagePixelType;
  typedef typename std::set<ImagePixelType>						PixelTypeSet;

  typedef typename itk::Vector<unsigned int,FEMDimension>							VectorTriIdType;
  typedef typename itk::MapContainer<unsigned int, PixelTypeSet>					SurfaceVertexToLabelMapContainerType;
  typedef typename itk::MapContainer<unsigned int, std::set<unsigned int>* > 		SurfaceVertexToNeighborVerticesMapContainerType;
  typedef typename itk::MapContainer<unsigned int, VectorType>						SurfaceVertexToNormalMapContainerType;
  typedef typename itk::MapContainer<unsigned int, VectorType >						SurfaceVertexToCorrespondenceMapContainerType;
  typedef typename itk::VectorContainer<unsigned int,VectorTriIdType>				SurfaceTrianglesContainerType;
  typedef typename itk::PointSet<PixelTypeSet>										TargetPointSetType;
  typedef typename itk::MapContainer<unsigned int, std::vector<unsigned int>* >		SurfaceVertexToNeighborTargetPointsMapContainerType;


  typedef typename itk::VectorContainer<unsigned int, ImagePixelType>				TetrahedraDataContainerType;

  typedef typename TargetPointSetType::PointsContainerIterator						TargetPointSetIteratorType;
  typedef typename TargetPointSetType::PointDataContainerIterator					TargetPointDataSetIteratorType;
  typedef typename Superclass::FEMObjectType    									FEMObjectType;

  /** Some convenient types */
  typedef typename Element::MatrixType                     							MatrixType;
  typedef typename FEMObjectType::LoadContainerType        							LoadContainerType;
  typedef typename FEMObjectType::NodeContainerType        							NodeContainerType;
  typedef typename FEMObjectType::LoadContainerIterator    							LoadContainerIterator;

  typedef LoadNoisyLandmark 														LoadType;

  /**
   * Number of iterations for the compressor.
   */
  itkSetMacro(NumIterations, unsigned int);
  itkGetMacro(NumIterations, unsigned int);

  itkSetMacro(ScaleFactorForEdges, double);
  itkGetMacro(ScaleFactorForEdges, double);

  itkSetMacro(ScaleFactorForVolumes, double);
  itkGetMacro(ScaleFactorForVolumes, double);

  itkSetMacro(ScaleFactorForSearchRadius, double);
  itkGetMacro(ScaleFactorForSearchRadius, double);

   /**
   * Control the balance between the Mesh deformation energy and the matching
   * energy. The range of values is from 0.0 to 1.0. When set to zero, the matching
   * energy is not considered. When set to 1.0 the solver will consider equally the
   * Mesh energy and the matching energy.
   */
  itkSetMacro(Lambda, double);
  itkGetMacro(Lambda, double);

  itkSetObjectMacro(SurfaceVertexToLabelMapContainer, SurfaceVertexToLabelMapContainerType);
  itkGetObjectMacro(SurfaceVertexToLabelMapContainer, SurfaceVertexToLabelMapContainerType);

  itkSetObjectMacro(SurfaceVertexToNeighborVerticesMapContainer, SurfaceVertexToNeighborVerticesMapContainerType);
  itkGetObjectMacro(SurfaceVertexToNeighborVerticesMapContainer, SurfaceVertexToNeighborVerticesMapContainerType);

  itkSetObjectMacro(SurfaceTrianglesContainer, SurfaceTrianglesContainerType);
  itkGetObjectMacro(SurfaceTrianglesContainer, SurfaceTrianglesContainerType);

  itkSetObjectMacro(TargetPointSet, TargetPointSetType);
  itkGetObjectMacro(TargetPointSet, TargetPointSetType);

  itkSetObjectMacro(TetrahedraDataContainer, TetrahedraDataContainerType);
  itkGetObjectMacro(TetrahedraDataContainer, TetrahedraDataContainerType);

protected:

  /**
   * Default constructor which sets the indices
   * for the matrix and vector storage.
   */
  CompressorSolver();
  ~CompressorSolver();

//  /** Method invoked by the pipeline in order to trigger the computation of
//   * the registration. */
  void GenerateData();

  void Initialization();

  void AddSurfaceVertexToItsOwnNeighborsSet();
  void RemoveSurfaceVertexFromItsOwnNeighborsSet();

  /**
   * Assemble the global mechanical stiffness matrix from the mesh contained in
   * the FEMObject
   */
  void AssembleMeshStiffnessMatrix();

  /**
   * Assemble element stiffness matrix, which will be used to assemble the
   * global stiffness matrix
   */
  virtual void AssembleElementMatrixWithID(const Element::Pointer & e, unsigned int matrixIndex);

  /**
   * Compute the tensor associated with the mesh surface vertices. The tensor is structural
   * weighted if a structural tensor point set is available
   */
  void ComputeLandmarksTensor();

  /**
   * Get scaling factor
   */
  float GetLandmarkTensorPonderation() const;

  /**
   * Simulate the landmark as a physical point and
   * assemble its contribution matrix
   */
  void AssembleLandmarkStiffnessMatrix();

  /**
   * Add global stiffness matrix with landmark stiffness matrix
   */
  void AssembleGlobalMatrixFromLandmarksAndMeshMatrices();

  /**
  * Assemble right side F vector based on the landmarks
  */
  void AssembleF();

  /**
  * Compute the surface vertices (landmarks) normals
  */
  void ComputeNormals();

  /**
   * Compute the Correspondences between source and target points
   */
  void ComputeCorrespondence();

  void UpdateLandmarks();

  /**
   * Compute the normal of a triangle
   */
  VectorType ComputeTriangleNormal(VectorTriIdType id);

  /**
   * Normalize all computed normals
   */
  void NormalizeNormals();


  void PrintForces(const char* fName);
  void PrintSolution(const char * fName);
  void UpdateFEMObjectCoordinates();
  void RestorePreviousFEMObjectCoordinates();
  void EvaluateCorrespondenceWithEdgeLength(VectorType& dVec,unsigned int vId,const double fac);
  void EvaluateCorrespondenceWithVolume(VectorType& dVec,VectorType& newPos,unsigned int vId,const double fac);
  void ComputeNeighborTargetPoints(const double fac);
  void PrintNeighborTargetVertices(std::vector<unsigned int>* pVector,typename TargetPointSetType::PointsContainer* pPointsContainer,std::string fName);

  bool IsInsideBox(typename TargetPointSetType::PointType pt,VectorType minCoord,VectorType maxCoord)
  {
	  for(unsigned int i=0; i < FEMDimension; i++)
	  {
		  if(pt[i] < minCoord[i] || pt[i] > maxCoord[i])
			  return false;
	  }
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

  bool PartiallySameLabelSets(PixelTypeSet set1,PixelTypeSet set2)
  {
	  if(set1.size() != 2 || set2.size() != 2)
		  /*return false;*/
		  itkExceptionMacro("set1 or set2 have size different than 2.");

	  ImagePixelType pix1,pix2,pix3,pix4;
	  typename PixelTypeSet::iterator it1 = set1.begin();
	  pix1 = *it1; ++it1;	pix2 = *it1;
	  typename PixelTypeSet::iterator it2 = set2.begin();
	  pix3 = *it2; ++it2;	pix4 = *it2;

	  if( pix1 == pix3  || pix1 == pix4 )
		  return true;
	  if( pix2 == pix3  || pix2 == pix4 )
		  return true;

	  return false;
  }

  bool CheckTetrahedraVolumes(itk::VectorContainer<unsigned long ,itk::fem::Element::Pointer> *ElementContainer)
  {
	  std::cout << "Checking tetrahedra volumes..." << std::endl;
	  const double ZERO6 = 0.000001;

	  if(ElementContainer == NULL)
		  return false;

	  int numTets = ElementContainer->Size();
	  itk::fem::Element::Pointer pTetr = NULL;

	  int count=0;

	  vnl_vector<double> pos0,pos1,pos2,pos3;
	  double vol=0;
	  for( int i=0; i < numTets; i++ )
	  {
		  pTetr = ElementContainer->at(i);
		  if(pTetr.IsNull())
		  {
			  std::cout << "Found NULL tetrahedron!" << std::endl;
			  exit(-1);
		  }

		  pos0 = pTetr->GetNode(0)->GetCoordinates();
		  pos1 = pTetr->GetNode(1)->GetCoordinates();
		  pos2 = pTetr->GetNode(2)->GetCoordinates();
		  pos3 = pTetr->GetNode(3)->GetCoordinates();

		  vol = CalculateTetrahedronVolume(pos0,pos1,pos2,pos3);
		  if(vol < ZERO6)
			  count++;
	  }

	  if(count == 0)
		  return true;
	  else
	  {
		  std::cout << count << " inverse tetrahedra found!" << std::endl;
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


private:

  CompressorSolver(const Self &);   //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  unsigned int m_NumIterations;

  /**
   * This type represents the index of the vector and matrix.
   */
  typedef unsigned int FEMIndexType;

  FEMIndexType 	m_ForceIndex;
  FEMIndexType 	m_SolutionIndex;
  FEMIndexType 	m_MeshStiffnessMatrixIndex;
  FEMIndexType 	m_LandmarkStiffnessMatrixIndex;
  FEMIndexType 	m_StiffnessMatrixIndex;

  float	   		m_Lambda;
  double 		m_ScaleFactorForEdges;
  double 		m_ScaleFactorForVolumes;
  double 		m_ScaleFactorForSearchRadius;

  LinearSystemWrapperItpack m_Itpack;

  typename SurfaceVertexToLabelMapContainerType::Pointer				m_SurfaceVertexToLabelMapContainer;
  typename SurfaceVertexToNormalMapContainerType::Pointer				m_SurfaceVertexToNormalMapContainer;
  typename SurfaceVertexToCorrespondenceMapContainerType::Pointer		m_SurfaceVertexToCorrespondenceMapContainer;
  typename SurfaceVertexToNeighborVerticesMapContainerType::Pointer 	m_SurfaceVertexToNeighborVerticesMapContainer;
  typename SurfaceTrianglesContainerType::Pointer						m_SurfaceTrianglesContainer;
  typename TargetPointSetType::Pointer									m_TargetPointSet;
  typename SurfaceVertexToNeighborTargetPointsMapContainerType::Pointer m_SurfaceVertexToNeighborTargetPointsMapContainer;
  ;
  // Helper container for printing the cell data
  typename TetrahedraDataContainerType::Pointer							m_TetrahedraDataContainer;


};

}

}  // end namespace itk::fem

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFEMCompressorSolver.hxx"
#endif

#endif /* ITKFEMCOMPRESSORSOLVER_H_ */
