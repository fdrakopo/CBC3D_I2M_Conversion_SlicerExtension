Extension: Image-To-Mesh Conversion (I2M) for Image Guided Therapy.

Description: This Slicer extension encapsulates two CLI modules: (1) Body Centric Cubic (BCC) Mesh Generation. This module generates a Body Centric Cubic (BCC) mesh from a labeled image. Initially the generated mesh is homogeneous, that means does not distinguish different tissues. Later the component specifies which tissue each tetrahedron belongs to. Each tissue is capable of automatically adjusting its resolution based on its geometric complexity and the predefined subdivision criterion. (2) Mesh Compression (MC). This module deforms an input tetrahedral mesh towards the boundaries of the input labeled image. Two point sets are extracted for the mesh deformation. The first (source point set) consists of the surface vertices of the input mesh. The second (target point set) consists of the surface edge points in the input labeled image. Then the input mesh is deformed by registering the source to the target point set using a Physics-Based Non-Rigid Registration method.

Contributors: Fotis Drakopoulos (CRTC), Yixun Liu (CRTC), Andrey Fedorov (CRTC), Nikos Chrisochoides (CRTC)
Center for Real-Time Computing (https://crtc.cs.odu.edu) &
Surgical Planning Lab, Brigham and Women's Hospital (http://www.spl.harvard.edu/)
Old Dominion Univeristy and The College of William and Mary, Virginia, US.
Contact:  npchris@gmail.com

Acknowledgements: This work was partially supported by NIH R44 OD018334-03A,
NSF grant No. CCF-1139864 and by the Richard T.Cheng Endowment.

Wiki Documentation: http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/ImageToMeshConversion

References: 
1. Tetrahedral mesh generation for medical imaging Fedorov A., 
Chrisochoides N., Kikinis R., Warfield S.
College of William and Mary; Computational Radiology Lab and Surgical 
Planning Lab, Brigham and Women's Hospital, Harvard Medical School

2. Mesh Deformation-based Multi-tissue Mesh Generation for Brain Images. 
Yixun Liu, Panagiotis Foteinos, Andrey Chernikov and Nikos Chrisochoides. 
In Engineering with Computers, Volume 28, pages 305-318, 2012.
