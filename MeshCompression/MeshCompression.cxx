#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"
#include "itkTimeProbesCollectorBase.h"

#include "MeshCompressionCLP.h"
#include "itkVTKUnstructuredGridReader.h"
#include "itkFEMCompressor.h"
#include "itkTetrahedralMeshWriter.h"


// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{


template <class T>
int DoIt( int argc, char * argv[], T )
{
	PARSE_ARGS;

	// Add a time probe
	itk::TimeProbesCollectorBase collector;

	const    unsigned int ImageDimension = 3;
	typedef    T InputPixelType;

	typedef itk::Image<InputPixelType,ImageDimension> InputImageType;
	typedef itk::ImageFileReader<InputImageType>  ReaderType;

	typename ReaderType::Pointer pImageReader = ReaderType::New();
	pImageReader->SetFileName( InputLabeledImageFileName.c_str() );
	try
	{
	    collector.Start( "Read input labeled image" );
		pImageReader->Update();
	    collector.Stop( "Read input labeled image" );
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error Reading Labeled Image: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	typename InputImageType::Pointer pInputImage = pImageReader->GetOutput();

	// Prints
	//std::cout << "Input Mesh : " << InputMeshFileName << std::endl;
	//std::cout << "Labeled Image : " << InputLabeledImageFileName << std::endl;
	//std::cout << "Trade Off : "  << Flexibility << std::endl;
	//std::cout << "Number of Iterations : "  << NumOfIterations << std::endl;
	//std::cout << "Non-Connectivity : " << NonConnectivity << std::endl;
	//std::cout << "Output Mesh : " << OutputMeshFileName << std::endl;
	//std::cout << "Output Source Points : " << OutputSourcePointsFileName << std::endl;
	//std::cout << "Output Target Points : " << OutputTargetPointsFileName << std::endl;
	//std::cout << "Output Surface Triangles : " << OutputSurfaceTrianglesFileName << std::endl;
	//std::cout << "" << std::endl;

	// Check parameters
    // (0.1 , 1]
    if(Flexibility < 0.1 || Flexibility > 1.0)
    {
        std::cerr << "Invalid value for Flexibility :  " << Flexibility
                  << "\nShould be between (0.1 , 1]" << std::endl;
        return EXIT_FAILURE;
    }
    // [1,10]
    if(NumOfIterations < 0 || NumOfIterations > 10)
    {
        std::cerr << "Invalid value for NumOfIterations :  " << NumOfIterations
                  << "\nShould be between [1,10]." << std::endl;
        return EXIT_FAILURE;
    }
    // [0,3]
	if(NonConnectivity != 0 && NonConnectivity != 1 && NonConnectivity != 2 && NonConnectivity != 3)
	{
	    std::cerr << "Invalid value for Non-Connectivity :  " << NonConnectivity
	              << "\nShould be 0 (vertex non-Connectivity) or 1 (edge non-Connectivity) or 2 (face non-Connectivity) or 3 (no non-Connectivity)" << std::endl;
	    return EXIT_FAILURE;
	}

	// Read the input mesh
	std::cout << "\nReading input mesh..." << std::endl;
	typedef float MeshPixelType;
	typedef itk::Mesh<MeshPixelType,ImageDimension> MeshType;
	typename itk::VTKUnstructuredGridReader<MeshType>::Pointer pMeshReader = itk::VTKUnstructuredGridReader<MeshType>::New();
	pMeshReader->SetFileName(InputMeshFileName);
	try
	{
	    collector.Start( "Read input mesh" );
		pMeshReader->Update();
	    collector.Stop( "Read input mesh" );
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "Error reading input mesh: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	MeshType::Pointer pInputMesh = pMeshReader->GetOutput();
    std::cout << "Number of Tetrahedra : " << pInputMesh->GetNumberOfCells() << std::endl;
    std::cout << "Number of Vertices : " << pInputMesh->GetNumberOfPoints() << std::endl;
	//std::cout << "Number of Cell Data: " << pInputMesh->GetCellData()->Size() << std::endl;

	// Call the FEMCompressor filter
	typedef itk::fem::FEMCompressor<MeshType,InputImageType,MeshType> FEMCompressorType;
	typename FEMCompressorType::Pointer pMeshCompressor = FEMCompressorType::New();
	pMeshCompressor->SetInput(pInputMesh);
	pMeshCompressor->SetImage(pInputImage);
	pMeshCompressor->SetLambda(Flexibility);
	pMeshCompressor->SetNumIterations(NumOfIterations);
	pMeshCompressor->SetNonConnectivityType(NonConnectivity);
	pMeshCompressor->SetOutputSourcePointsFileName(OutputSourcePointsFileName);
	pMeshCompressor->SetOutputTargetPointsFileName(OutputTargetPointsFileName);
	pMeshCompressor->SetOutputSurfaceTrianglesFileName(OutputSurfaceTrianglesFileName);
	try
	{
	    collector.Start( "Compress mesh" );
		pMeshCompressor->Update();
	    collector.Stop( "Compress mesh" );
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "Error in Mesh Compressor: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	// Write the output mesh to a vtk file
	if(!OutputMeshFileName.empty())
	{
		itk::TetrahedralMeshWriter<MeshType>::Pointer pMeshWriter = itk::TetrahedralMeshWriter<MeshType>::New();
		pMeshWriter->SetFileName(OutputMeshFileName.c_str());
		pMeshWriter->SetMesh(pMeshCompressor->GetOutput());
		try
		{
		    collector.Start( "Write mesh" );
			pMeshWriter->Update();
		    collector.Stop( "Write mesh" );
		}
		catch(itk::ExceptionObject &err)
		{
			std::cerr << "Failed to write output mesh : " << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}

	collector.Report();
	return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
	//InputLabeledImageFileName = "/home/fdrakopo/DEVELOPMENT/ImageToMeshConversion_Slicer_MultiTissue/I2M/MeshCompression/Data/Baseline/nidus.mha";
    itk::GetImageType(InputLabeledImageFileName, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
