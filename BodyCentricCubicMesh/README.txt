BCC_I2M_Tool - Body-Centered Cubic (BCC) lattice-based Image-to-Mesh Conversion Tool

Contributors: Fotis Drakopoulos, Yixun Liu, Andrey Fedorov, Nikos Chrisochoides
Center for Real-Time Computing (https://crtc.cs.odu.edu) &
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

The method generates a Body Centric Cubic (BCC) mesh from a labeled
image. Initially the generated mesh is homogeneous, that means does not
distinguish different tissues. Later the module specifies which tissue
each tetrahedron belongs to. Each tissue is capable of automatically
adjusting its resolution based on its geometric complexity and the
predefined subdivision criterion. The current version of the module
supports a single-tissue input labeled image.


Supported Platforms 
Linux

Building BCC_I2M_Tool
This tool was built as a CLI (Command Line Interface) 3D Slicer module. 
A stand alone ITK version of the code is also available. For more information please contact Nikos Chrisochoides (npchris@gmail.com).


USAGE: 

   lib/Slicer-4.3/cli-modules/BodyCentricCubicMesh  [--returnparameterfile
                                        <std::string>]
                                        [--processinformationaddress
                                        <std::string>] [--xml] [--echo]
                                        [--outputMesh <std::string>]
                                        [--inputLabeledImage <std::string>]
                                        [--resampleResolution <int>]
                                        [--fidelity <double>] [--size
                                        <double>] [--] [--version] [-h]


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

   --outputMesh <std::string>
     The output tetrahedral mesh in vtk file format. (default:
     outputMesh.vtk)

   --inputLabeledImage <std::string>
     The input labeled image.

   --resampleResolution <int>
     Controls the image resampling. Suggested value is 0 (no-resampling)
     for large tissues and 1 or 2 for smaller tissues. Value should be 0 or
     1 or 2. (default: 0)

   --fidelity <double>
     Controls the subdivision of the elements that belong to more than one
     image labels. The larger the value the better the conformity of the
     mesh on the image label interfaces. Value should be between [0.1,1].
     (default: 0.8)

   --size <double>
     Controls the size of the mesh. The smaller the value the larger the
     mesh. Value should be between [1,20]. (default: 8)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Description: Generates a Body Centric Cubic (BCC) mesh from a labeled
   image. Initially the generated mesh is homogeneous, that means does not
   distinguish different tissues. Later the module specifies which tissue
   each tetrahedron belongs to. Each tissue is capable of automatically
   adjusting its resolution based on its geometric complexity and the
   predefined subdivision criterion. The current version of the module
   supports a single-tissue input labeled image. A comprehensive
   description of the method can be found at : 1.Tetrahedral mesh
   generation for medical imaging Fedorov A., Chrisochoides N., Kikinis R.,
   Warfield S. College of William and Mary; Computational Radiology Lab and
   Surgical Planning Lab, Brigham and Women's Hospital, Harvard Medical
   School, and 2. Mesh Deformation-based Multi-tissue Mesh Generation for
   Brain Images, Yixun Liu, Panagiotis Foteinos, Andrey Chernikov and Nikos
   Chrisochoides, Engineering with Computers, Volume 28, pages 305-318,
   2012.

   Author(s): Fotis Drakopoulos (CRTC), Yixun Liu (CRTC), Andrey Fedorov
   (CRTC), Nikos Chrisochoides (CRTC)

   Acknowledgements: This work was partially supported by NIH R44
   OD018334-03A, grant CCF-1139864 and by the Richard T.Cheng Endowment.

