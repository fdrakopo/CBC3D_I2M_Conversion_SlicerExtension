cmake_minimum_required(VERSION 2.8.9)

project(ImageToMeshConversion)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/ImageToMeshConversion")
set(EXTENSION_CATEGORY "Mesh Generation")
set(EXTENSION_CONTRIBUTORS "Fotis Drakopoulos (CRTC), Yixun Liu (CRTC), Andrey Fedorov (CRTC), Nikos Chrisochoides (CRTC)")
set(EXTENSION_DESCRIPTION "Image-To-Mesh Conversion (I2M) for image Guided Therapy. The extension consists of two modules: (1) Body Centric Cubic (BCC) Mesh Generation.This module generates a Body Centric Cubic (BCC) mesh from a labeled image. Initially the generated mesh is homogeneous, that means does not distinguish different tissues. Later the component specifies which tissue each tetrahedron belongs to. Each tissue is capable of automatically adjusting its resolution based on its geometric complexity and the predefined subdivision criterion. (2) Mesh Compression (MC). This module deforms an input tetrahedral mesh towards the boundaries of the input labeled image. Two point sets are extracted for the mesh deformation. The first (source point set) consists of the exterior and interior surface vertices of the input mesh. The second (target point set) consists of the exterior and interior surface edge points in the input labeled image. Then the input mesh is deformed by registering the source to the target point set using a Physics-Based Non-Rigid Registration method.")
set(EXTENSION_ICONURL "http://www.slicer.org/slicerWiki/images/e/ef/CRTC_logo.jpg")
set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/e/ea/Brain_before.jpg http://www.slicer.org/slicerWiki/images/b/b4/Brain_after.jpg http://www.slicer.org/slicerWiki/images/f/f3/Brain_after_section.jpg http://www.slicer.org/slicerWiki/images/5/54/Ventricles_before.jpg http://www.slicer.org/slicerWiki/images/3/31/Ventricles_after.jpg http://www.slicer.org/slicerWiki/images/7/7c/Ventricles_after_section.jpg")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(BodyCentricCubicMesh)
add_subdirectory(MeshCompression)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
