/*=========================================================================


=========================================================================*/
#ifndef __itkTetrahedralMeshWriter_h
#define __itkTetrahedralMeshWriter_h

#include "itkMesh.h"
#include "itkCellInterface.h"
#include "itkTetrahedronCell.h"
#include "itkMeshSource.h"

namespace itk
{

/** \class itkTetrahedralMeshWriter
 * 
 * 
 * \par
 * This class writes an ITK mesh into a vtk file
 * The mesh should be tetrahedral. 
 *
 * \par PARAMETERS
 * 
 * \par REFERENCE
 *
 * \par INPUT
 * ITK tetrahedral mesh
 *
 * */
  
template <class TInputMesh>
class ITK_EXPORT TetrahedralMeshWriter: public MeshSource<TInputMesh>
{
public:
  typedef TetrahedralMeshWriter Self;
  typedef MeshSource<TInputMesh> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(TetrahedralMeshWriter, MeshSource);

  typedef TInputMesh MeshType;
  typedef typename MeshType::MeshTraits   		MeshTraits;
  typedef typename MeshType::CellTraits   		CellTraits;
  typedef typename MeshType::PointType    		PointType;
  typedef typename MeshType::CellType     		CellType;
  typedef typename MeshTraits::PixelType  		PixelType;
  typedef typename CellType::PointIdIterator  	PointIdIterator;
  
  typedef typename MeshType::PointsContainer 		PointsContainer;
  typedef typename MeshType::PointsContainerPointer PointsContainerPointer;
  typedef typename MeshType::CellsContainer         CellsContainer;
  typedef typename MeshType::CellDataContainer      CellDataContainer;


  typedef CellInterface<PixelType, CellTraits> 	TCellInterface;
  typedef TetrahedronCell<TCellInterface> 		TetrahedronCellType;

  itkSetMacro(FileName, std::string);
  void SetMesh(MeshType* pMesh)
  {
	  this->m_InputMesh = pMesh;
  }

protected:
  TetrahedralMeshWriter();
  ~TetrahedralMeshWriter();

private:
  typename MeshType::Pointer m_InputMesh;
  std::string m_FileName;

  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData();
  void TestFileExistanceAndReadability();
  void WriteVTKMesh();
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTetrahedralMeshWriter.cxx"
#endif 

#endif // __itkTetrahedralMeshWriter_h
