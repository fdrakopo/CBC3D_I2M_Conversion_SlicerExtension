#ifndef __itkTetrahedralMeshWriter_txx
#define __itkTetrahedralMeshWriter_txx

#include "itkTetrahedralMeshWriter.h"
#include <itksys/SystemTools.hxx>

#include "vtkUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkFloatArray.h"
#include "vtkCellData.h"

#include "vtkNew.h"


namespace itk
{

template<class TInputMesh>
TetrahedralMeshWriter<TInputMesh>::TetrahedralMeshWriter()
{
	this->m_FileName = "";
	this->m_InputMesh = NULL;
}

template<class TInputMesh>
TetrahedralMeshWriter<TInputMesh>::~TetrahedralMeshWriter()
{
}

template<class TInputMesh>
void TetrahedralMeshWriter<TInputMesh>::PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf(os, indent);
	os << indent << "m_FileName: " << m_FileName << "\n";
}


template<class TInputMesh>
void TetrahedralMeshWriter<TInputMesh>::TestFileExistanceAndReadability()
{
	// Test if the file exists.
	if( ! itksys::SystemTools::FileExists( m_FileName.c_str() ) )
	{
		itkExceptionMacro(<<"File does not exist");
	}

	// Test if the file can be open for reading access.
	std::ifstream readTester;
	readTester.open( this->m_FileName.c_str() );
	if( readTester.fail() )
	{
		readTester.close();
		itkExceptionMacro(<<"File cannot be open");
		return;
	}
	readTester.close();
}

template<class TInputMesh>
void TetrahedralMeshWriter<TInputMesh>::GenerateData()
{
	if(this->m_FileName == "")
		itkExceptionMacro(<<"Missing input file name");

	if(this->m_InputMesh.IsNull())
		itkExceptionMacro(<<"The input mesh is null");

	// Find out the type of data saved
	std::string suffix = m_FileName.substr(m_FileName.rfind(".")+1,m_FileName.size());
	this->m_InputMesh->SetBufferedRegion(this->m_InputMesh->GetRequestedRegion());

	if(suffix == "vtk")
		WriteVTKMesh();
	else
		itkExceptionMacro(<<"File type not supported");

	this->SetNthOutput(0, this->m_InputMesh);
}

template<class TInputMesh>
void TetrahedralMeshWriter<TInputMesh>::WriteVTKMesh()
{
	std::cout << "Writing output mesh at : " << this->m_FileName << std::endl;

	if(this->m_InputMesh.IsNull())
		itkExceptionMacro(<<"The input mesh is null");

	// Get the points in the mesh
	PointsContainerPointer points = this->m_InputMesh->GetPoints();
	int numPoints = points->Size();
	if(numPoints == 0)
		itkExceptionMacro(<<"The number of mesh points is zero.");

	vtkUnstructuredGrid* vgrid = vtkUnstructuredGrid::New();
	vtkPoints* vpoints = vtkPoints::New();
	vpoints->SetNumberOfPoints(numPoints);

	// Set the points
	typename PointsContainer::ConstIterator pointsIter = points->Begin();
	unsigned int id = 0;
	PointType pt;

	while(pointsIter != points->End())
	{
		id = pointsIter->Index();
		for(unsigned int i=0; i<3; i++)
			pt[i] = pointsIter->Value()[i];

		vpoints->SetPoint(id, pt[0], pt[1], pt[2]);
		++pointsIter;
	}
	vgrid->SetPoints(vpoints);
	vpoints->Delete();

	// Set the cells and the data
	vtkIdList* ptIds = vtkIdList::New();
	vtkFloatArray* labelArray = vtkFloatArray::New();
	labelArray->SetNumberOfComponents(1);
	labelArray->SetNumberOfTuples(this->m_InputMesh->GetNumberOfCells());
	labelArray->SetName("Label");

	if(!this->m_InputMesh->GetCells())
		itkExceptionMacro("The Cells container is null.");
	if(!this->m_InputMesh->GetCellData())
		itkExceptionMacro("The Cell Data container is null.");
	if(this->m_InputMesh->GetCellData()->Size() != this->m_InputMesh->GetCells()->Size())
		itkExceptionMacro("The nunber of Cells data is different than the number of Cells.");

	typename CellsContainer::ConstIterator cellIter = this->m_InputMesh->GetCells()->Begin();
	typename CellsContainer::ConstIterator cellEndIter = this->m_InputMesh->GetCells()->End();
	typename CellDataContainer::ConstIterator cellDataIter = this->m_InputMesh->GetCellData()->Begin();

	while( cellIter != cellEndIter )
	{
		CellType * pCell = cellIter.Value();
		if(!pCell)
			itkExceptionMacro(<<"The cell is null.");

		if(pCell->GetType() != CellType::TETRAHEDRON_CELL )
			itkExceptionMacro(<<"Invalid Cell Type. Type should be TETRAHEDRON_CELL.");

		TetrahedronCellType * pTetrCell = dynamic_cast<TetrahedronCellType*>(pCell);

		PointIdIterator pointIdIter = pTetrCell->PointIdsBegin();
		PointIdIterator pointIdEnd = pTetrCell->PointIdsEnd();

		unsigned int i=0;
		while( pointIdIter != pointIdEnd )
		{
			ptIds->InsertId(i,*pointIdIter);
			i++;
			++pointIdIter;
		}

		// Insert the cell to Unstructured Grid (10 = VTK_TETRA)
		vgrid->InsertNextCell(10,ptIds);

		// The label of the cell
		labelArray->InsertTuple1(cellDataIter.Index(), cellDataIter.Value());

		++cellIter;
		++cellDataIter;
	}

	ptIds->Delete();
	vgrid->GetCellData()->AddArray(labelArray);
	labelArray->Delete();

	vtkNew<vtkUnstructuredGridWriter> vtk_output_writer;
	vtk_output_writer->SetFileName(m_FileName.c_str());
	vtk_output_writer->SetInputData(vgrid);
	vtk_output_writer->Update();
	vgrid->Delete();
}


}

#endif // __itkTetrahedralMeshWriter_txx












