#include "itkPluginUtilities.h"

#include "itkTimeProbesCollectorBase.h"

#include "BodyCentricCubicMeshCLP.h"
#include "itkBinaryMaskTo3DAdaptiveMeshFilter.h"
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

	pImageReader->SetFileName( LabeledImageFileName.c_str() );
	try
	{
	    collector.Start( "Read input image" );
		pImageReader->Update();
	    collector.Stop( "Read input image" );
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error Reading Labeled Image: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	// Prints
	std::cout << "" << std::endl;
	std::cout << "Input image : " << LabeledImageFileName << std::endl;
	std::cout << "Mesh Size : " << MeshSize << std::endl;
	std::cout << "Fidelity : " << Fidelity << std::endl;
	std::cout << "Resample Resolution : "  << ResampleResolution << std::endl;
	std::cout << "Output Mesh File Name : "  << OutputMeshFileName << std::endl;
	std::cout << "" << std::endl;

	// Check parameters
	if(MeshSize < 1 || MeshSize > 20)
	{
	    std::cerr << "Invalid value for Mesh Size :  " << MeshSize
	              << "\nShould be between [1,20]." << std::endl;
	    return EXIT_FAILURE;
	}
	if(Fidelity < 0.1 || Fidelity > 1.0)
	{
	    std::cerr << "Invalid value for Fidelity :  " << Fidelity
	              << "\nShould be between [0.1, 1]." << std::endl;
	    return EXIT_FAILURE;
	}
	if(ResampleResolution != 0 && ResampleResolution != 1 && ResampleResolution != 2)
	{
	    std::cerr << "Invalid value for ResampleResolution :  " << ResampleResolution
	              << "\nShould be 0 or 1 or 2." << std::endl;
	    return EXIT_FAILURE;
	}

	typename InputImageType::Pointer pInputImage = pImageReader->GetOutput();
	typename InputImageType::SizeType input_size = pInputImage->GetLargestPossibleRegion().GetSize();
	typename InputImageType::SpacingType input_spacing = pInputImage->GetSpacing();
	std::cout<<"Input Image: " << std::endl;
	std::cout<<"Size : " << input_size << std::endl;
	std::cout<<"Spacing : " << input_spacing << std::endl;
	std::cout<<"Origin : " << pInputImage->GetOrigin()<< std::endl;
	std::cout<<"Direction : " << pInputImage->GetDirection()<< std::endl;
	std::cout << "" << std::endl;

	// This variable defines the spacing of the initial BCC lattice. By default, it is equal to 10,
	// which means that the spacing will be 1/10th of the smallest dimension of the image
	int BCCSpacing = fmin(fmin(input_size[0]*input_spacing[0],input_size[1]*input_spacing[1]), input_size[2]*input_spacing[2]) / MeshSize;

	// The mesh filter
	typedef float MeshPixelType;
	typedef itk::Mesh<MeshPixelType,ImageDimension> MeshType;
	typedef itk::BinaryMaskTo3DAdaptiveMeshFilter<InputImageType,MeshType> MeshFilterType;
	typename MeshFilterType::Pointer pMeshFilter = MeshFilterType::New();

	pMeshFilter->SetInput(pInputImage);
	pMeshFilter->SetSubdivInclusion(true);
	pMeshFilter->SetSubdivCurvature(false);
	pMeshFilter->SetBCCSpacing(BCCSpacing);
	pMeshFilter->SetFidelity(Fidelity);
	pMeshFilter->SetResampleResolution(ResampleResolution);
	try
	{
	    collector.Start( "Generate mesh" );
		pMeshFilter->Update();
	    collector.Stop( "Generate mesh" );
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "Failed to update the mesh filter : " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	// Write the generated tetrahedral mesh to a vtk file
	if(!OutputMeshFileName.empty())
	{
		itk::TetrahedralMeshWriter<MeshType>::Pointer pMeshWriter = itk::TetrahedralMeshWriter<MeshType>::New();
		pMeshWriter->SetFileName(OutputMeshFileName.c_str());
		pMeshWriter->SetMesh(pMeshFilter->GetOutput());
		try
		{
		    collector.Start( "Write output mesh" );
			pMeshWriter->Update();
		    collector.Stop( "Write output mesh" );
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
	  itk::GetImageType(LabeledImageFileName, pixelType, componentType);

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
