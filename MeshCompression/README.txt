BCC_I2M_Tool - Mesh Compression (MC) Image-to-Mesh Conversion Tool

Contributors: Fotis Drakopoulos, Yixun Liu, Andrey Fedorov, Nikos Chrisochoides
Center for Real-Time Computing (crtc.wm.edu) &
Surgical Planning Lab, Brigham and Women's Hospital (http://www.spl.harvard.edu/)
Old Dominion Univeristy and The College of William and Mary, Virginia, US.
Contact:  npchris@gmail.com

Project Page: TBA

Acknowledgements: This work was partially supported by NIH R44 OD018334-03A, 
NSF grant No. CCF-1139864 and by the Richard T.Cheng Endowment.
 

The BCC_I2M_Tool is based on:
1. PhD thesis of Neil Molino (Stanford University, 2004
2. PhD Thesis of Andrey Fedorov (The College of William and Mary, 2009)
3. PhD Thesis of Yixun Liu (The College of William and Mary, 2011) 
4. PhD Thesis of Fotis Drakopoulos (Old Dominion University, in progress)
The PhD thesis 2 to 4 are supervised by Nikos Chrisochoides

The method for the BCC_I2M_Tool is described in:

1. Tetrahedral mesh generation for medical imaging Fedorov A., 
Chrisochoides N., Kikinis R., Warfield S.
College of William and Mary; Computational Radiology Lab and Surgical 
Planning Lab, Brigham and Women's Hospital, Harvard Medical School

2. Mesh Deformation-based Multi-tissue Mesh Generation for Brain Images. 
Yixun Liu, Panagiotis Foteinos, Andrey Chernikov and Nikos Chrisochoides. 
In Engineering with Computers, Volume 28, pages 305-318, 2012.

The Mesh Compression module deforms the input tetrahedral mesh towards the boundaries of the input labeled image. 
Two point sets are extracted for the mesh deformation. The first (source point set) consists of 
the surface vertices in the input mesh. The second (target point set) consists of the surface 
edge points in the input labeled image. After the extraction of the two point sets, 
the input mesh is deformed by registering the source to the target point set using a 
physics-based non-rigid registration method. The current version of the module supports a single-tissue input labeled image and mesh.

Supported Platforms 
Linux, Windows

Building BCC_I2M_Tool
This tool was built as a CLI (Command Line Interface) 3D Slicer module. 
The 3D Slicer tool needs to be downloaded, compiled and installed from source following the instructions here: http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Developers/Build_Instructions.
A stand alone ITK version of the code is also available. For more information please contact Nikos Chrisochoides (npchris@gmail.com).


USAGE: 

   lib/Slicer-4.3/cli-modules/MeshCompression  [--returnparameterfile
                                        <std::string>]
                                        [--processinformationaddress
                                        <std::string>] [--xml] [--echo]
                                        [--outputSurfaceTriangles
                                        <std::string>]
                                        [--outputTargetPoints
                                        <std::string>]
                                        [--outputSourcePoints
                                        <std::string>] [--outputMesh
                                        <std::string>] [--inputMesh
                                        <std::string>] [--inputLabeledImage
                                        <std::string>] [--nonConnectivity
                                        <int>] [--numIterations <int>]
                                        [--flexibility <double>] [--]
                                        [--version] [-h]


Where: 

   --returnparameterfile <std::string>
     Filename in which to write simple return parameters (int, float,
     int-vector, etc.) as opposed to bulk return parameters (image,
     geometry, transform, measurement, table).

   --processinformationaddress <std::string>
     Address of a structure to store process information (progress, abort,
     etc.). (default: 0)

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   --outputSurfaceTriangles <std::string>
     Extracted surface triangles of the input mesh in vtk file format
     (optional).

   --outputTargetPoints <std::string>
     Extracted target points in vtk file format (optional). The target
     points consist of the surface edge points and the interface edge
     points in the input labeled image.

   --outputSourcePoints <std::string>
     Extracted source points in vtk file format (optional). The source
     points consist of the surface vertices and the interface vertices
     (vertices that belong to more than one image label) in the input
     mesh.

   --outputMesh <std::string>
     Output deformed mesh in vtk file format (optional). (default:
     deformedMesh.vtk)

   --inputMesh <std::string>
     Input tetrahedral mesh in vtk file format.

   --inputLabeledImage <std::string>
     Input labeled image.

   --nonConnectivity <int>
     Controls the density of the extracted target points via a prohibited
     connectivity between adjacent blocks of voxels. The larger the value
     the denser the target point set. 0 is vertex non-connectivity (26
     neighbors), 1 is edge non-connectivity (18 neighbors), 2 is face
     non-connectivity (6 neighbors) and 3 is no non-connectivity. Value
     should be 0 or 1 or 2 or 3. (default: 0)

   --numIterations <int>
     Controls the number of iterations for the compressor. The larger the
     value the closer the boundaries of the input mesh to the boundaries of
     the input image. A large number of iterations might cause distorted
     poor-quality elements. Value should be between [1,10]. (default: 3)

   --flexibility <double>
     Controls the trade-off between the regularization energy defined by
     the stress energy of a linear biomechanical elastic model, and the
     similarity energy. Value should be between [0.1,1]. (default: 1)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Description: Deforms the input tetrahedral mesh towards the boundaries
   of the input labeled image. Two point sets are extracted for the mesh
   deformation. The first (source point set) consists of the surface
   vertices in the input mesh. The second (target point set) consists of
   the surface edge points in the input labeled image. After the extraction
   of the two point sets, the input mesh is deformed by registering the
   source to the target point set using a physics-based non-rigid
   registration method. The current version of the module supports a
   single-tissue input labeled image and mesh. A comprehensive description
   of the method can be found in article : Mesh deformation-based
   multi-tissue mesh generation for brain images, Engineering with
   Computers (2012), Volume 28, Issue 4, pp 305-318.

   Author(s): Fotis Drakopoulos (CRTC), Yixun Liu (CRTC), Andrey Fedorov
   (SPL B&W Harvard), Nikos Chrisochoides (CRTC)

   Acknowledgements: This work was partially supported by NIH R44
   OD018334-03A, grant CCF-1139864 and by the Richard T.Cheng Endowment.

