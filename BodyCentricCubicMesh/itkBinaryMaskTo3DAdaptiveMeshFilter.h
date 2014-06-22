#ifndef __itkBinaryMaskTo3DAdaptiveMeshFilter_h
#define __itkBinaryMaskTo3DAdaptiveMeshFilter_h

#include "itkImageToMeshFilter.h"
#include "itkTetrahedronCell.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"


extern "C"{
float exactinit();
double orient3d(double*, double*, double*, double*);
}

namespace itk
{

/** \class itkBinaryMaskTo3DAdaptiveMeshFilter
 * 
 * 
 * \par
 * This class tries to construct a 3D tetrahedral volume mesh based on a binary mask.
 *
 * \par
 * We start from a Body-Centered Cubic (BCC) lattice, with the initial BCC
 * spacing defined (TODO: use Medial Axis to detect the smallest feature size and
 * chose the spacing appropriately). Next selected tetrahedra are subdivided
 * (based on user-specified criteria). Neighbouring tetrahedra that have to be
 * subdivided to enforce the conformant mesh. The process is repeated up to the
 * user-specified number of resolutions.
 *
 * \par
 * After the maximum resolution has been reached, the mesh is modified to
 * ensure that the boundary is manifold, and to enforce some properties of the
 * surface, which proved to be useful to achieve good quality when compressing
 * the surface to the boundary.
 *
 * \par
 * The algorithm is driven by the distance function computed on the input binary
 * mask. The current subdivision criteria include the proximity of the surface
 * and the surface curvature both computed from the distance image. The
 * surface compression code is separate at the moment.
 * 
 * \par PARAMETERS
 * 
 *
 * \par REFERENCE
 *
 * Molino, N., Bridson, R., Teran, J., and Fedkiw, R. "A crystalline, red green 
 * strategy for meshing highly deformable objects with tetrahedra". 12th International 
 * Meshing Roundtable, Sandia National Laboratories, pp.103-114, Sept. 2003.
 * 
 * \par INPUT
 * The input should be a 3D binary image.
 *
 *  */
  
template <class TInputImage, class TOutputMesh>
class ITK_EXPORT BinaryMaskTo3DAdaptiveMeshFilter : public ImageToMeshFilter<TInputImage,TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef BinaryMaskTo3DAdaptiveMeshFilter         Self;
  typedef ImageToMeshFilter<TInputImage,TOutputMesh>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(BinaryMaskTo3DAdaptiveMeshFilter, ImageToMeshFilter);

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputImage InputImageType;
  typedef TOutputMesh OutputMeshType;
  typedef typename OutputMeshType::MeshTraits   OMeshTraits;
  typedef typename OutputMeshType::CellTraits   OCellTraits;
  typedef typename OutputMeshType::PointType    OPointType;
  typedef typename OutputMeshType::CellType     OCellType;
  typedef typename OMeshTraits::PixelType       OPixelType;

  /** Define tetrahedral cell type */
  typedef CellInterface<OPixelType, OCellTraits> TCellInterface;
  typedef TetrahedronCell<TCellInterface> TetCell;
  typedef typename TetCell::CellAutoPointer TetCellAutoPointer;

  /** Type definition of the subdivision test function */
  typedef bool (SubdivisionTestFunctionType)(double* v0, double* v1, double* v2, double* v3, Self*);
  typedef SubdivisionTestFunctionType* SubdivisionTestFunctionPointer;
  
  itkSetMacro(SubdivInclusion, bool);
  itkSetMacro(SubdivCurvature, bool);
  itkSetMacro(KeepOutside, bool);
  itkSetMacro(InputImagePrefix, std::string);
  itkSetMacro(TmpDirectory, std::string);
  itkGetMacro(NumberOfPoints, unsigned);
  itkGetMacro(NumberOfTets, unsigned);

  // Defines the spacing of the initial BCC lattice.
  // By default, it is equal to 10, which means that the spacing will be 1/10th of the
  // smallest dimension of the image
  itkSetMacro(BCCSpacing, unsigned);
  itkSetMacro(Fidelity, double);
  itkSetMacro(ResampleResolution, double);

protected:
  struct RGMVertex;
  struct RGMEdge;
  struct RGMTetra;
  typedef RGMVertex* RGMVertex_ptr;
  typedef RGMEdge*    RGMEdge_ptr;
  typedef RGMTetra*   RGMTetra_ptr;
  
  struct RGMVertex{
    double coords[3];
  };

  struct RGMEdge{
    RGMVertex_ptr             nodes[2];
    RGMEdge_ptr               children[2];  // cannot allocate two children at a time, because 
                                            // we do not keep all the levels
                                            // and children can be deallocated
                                            // arbitrarily
    RGMVertex_ptr             midv;
    std::list<RGMTetra_ptr>   neihood;
  };

  struct RGMTetra{            // 36 bytes
    RGMEdge_ptr               edges[6];
    RGMTetra_ptr              parent;
    RGMTetra_ptr              greens;       // green children; can be allocated as a consecutive region
    unsigned char             subdiv;
    unsigned char             edges_split;
    unsigned char             edges_orient;
    unsigned char             level;        
    unsigned char             num_greens;
	int						  label;		//for multi-tissues mesh generation
	RGMTetra_ptr s[4]; //for multi-tissues mesn generation. 4 possible shared Tetra. NULL denotes no share tetra.
  };

  	typedef struct{
		RGMTetra_ptr ptr0, ptr1;//two RGMTetra_ptrs surrounding the face.
	}FaceInfo;

  
  struct ltVertexPair{
    bool operator()(const std::pair<RGMVertex_ptr,RGMVertex_ptr> key1,
                    const std::pair<RGMVertex_ptr,RGMVertex_ptr> key2) const{
      if(key1.first==key2.first)
        return key1.second<key2.second;
      return key1.first<key2.first;
    }
  };
 
  
  BinaryMaskTo3DAdaptiveMeshFilter();
  ~BinaryMaskTo3DAdaptiveMeshFilter();
    
  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData();
  
private:
  class TetraRGFace{
  public:
    TetraRGFace(){};
    TetraRGFace(RGMVertex_ptr v0, RGMVertex_ptr v1, RGMVertex_ptr v2, RGMVertex_ptr v3){
      fourth = v3;
      if(v0<v1)
        if(v1<v2){
          nodes[0] = v0; nodes[1] = v1; nodes[2] = v2;
        } else {
          if(v2<v0){
            nodes[0] = v2; nodes[1] = v0; nodes[2] = v1;
          } else {
            nodes[0] = v0; nodes[1] = v2; nodes[2] = v1;
          }
        }
      else
        if(v0<v2){
          nodes[0] = v1; nodes[1] = v0; nodes[2] = v2;
        } else {
          if(v1<v2){
            nodes[0] = v1; nodes[1] = v2; nodes[2] = v0;
          } else {
            nodes[0] = v2; nodes[1] = v1; nodes[2] = v0;
          }
        }
    }
    
    ~TetraRGFace(){};
    
    bool operator==(const TetraRGFace &f) const{
      return ((f.nodes[0]==this->nodes[0]) && (f.nodes[1]==this->nodes[1]) &&
        (f.nodes[2]==this->nodes[2]));
    }
    
    bool operator<(const TetraRGFace &face1) const{
      if(this->nodes[0]<face1.nodes[0]) 
        return true;
      if(face1.nodes[0]<this->nodes[0])
        return false;
      if(this->nodes[1]<face1.nodes[1])
        return true;
      if(face1.nodes[1]<this->nodes[1])
        return false;
      return this->nodes[2]<face1.nodes[2];
    }

    RGMVertex_ptr nodes[3];
    RGMVertex_ptr fourth;
  };


  BinaryMaskTo3DAdaptiveMeshFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typedef float InternalPixelType;
  typedef Image<InternalPixelType,3> InternalImageType;
  typedef InternalImageType::PointType InternalImagePointType;
  typedef CastImageFilter<InputImageType,InternalImageType> CastFilterType;
  typedef NearestNeighborInterpolateImageFunction<InternalImageType,double> InterpolatorType;
  typedef ResampleImageFilter<InternalImageType,InternalImageType> ResamplerType;
  typedef IdentityTransform<double,3> IdentityTransformType;
  typedef ImageFileReader<InternalImageType> InternalImageReaderType;
  typedef ImageFileWriter<InternalImageType> InternalImageWriterType;

  // Lookup tables
  signed char m_ThirdFaceEdgeLT[36];
  unsigned char m_FaceToEdgesLT[12];
  unsigned char m_IncidentEdgesLT[24];
 
  double m_dimX, m_dimY, m_dimZ;
  double m_BCCSpacing;
  
  typedef bool HelperPixelType;
  typedef Image<HelperPixelType,3> HelperImageType;
  typename HelperImageType::Pointer m_HelperImage;

  std::set<int> m_labelSet;//the materials or labels needed to be meshed
  std::set<int> m_finalizedLabelSet;//for multi-tissues mesh generator
  std::map<int, int> m_labelToCount;//label to count map for the real region. Note that mesh region is different from real region

  std::list<RGMVertex_ptr> m_Vertices;  // these two lists may not be required.
  std::list<RGMEdge_ptr> m_Edges;       // I keep them for memory deallocation.
  std::list<RGMTetra_ptr> m_Tetras;     // (quick&dirty)
  std::list<RGMTetra_ptr> m_PendingTetras;
  std::map<std::pair<RGMVertex_ptr,RGMVertex_ptr>,RGMEdge_ptr,ltVertexPair> m_TmpEdgeMap;
  std::list<SubdivisionTestFunctionPointer> m_SubdivisionCriteria;
  
  typename InterpolatorType::Pointer m_Interpolator;
  typename InputImageType::Pointer m_InputImage;
  InternalImageType::Pointer m_ReadyInputImage;
  InternalImagePointType m_InputOrigin, m_OriginalInputOrigin;
  
  // Subdivision criteria
  bool m_SubdivInclusion;
  bool m_SubdivCurvature;
  
  bool m_KeepOutside;
  unsigned m_CurrentResol;
  std::string m_InputImagePrefix;
  std::string m_TmpDirectory;

  double m_Fidelity;
  int m_ResampleResolution;
  unsigned int m_MaxAdaptiveResolution;
  unsigned int m_MaxResolution;

  unsigned long m_NumberOfPoints;
  unsigned long m_NumberOfTets;

  bool Initialize();
  double FindBCCSpacing();
  void CreateMesh();
  void CreateBCC();
  bool TetraUnderrefined(RGMTetra_ptr);
  unsigned char GetTetraEdgeConf(RGMTetra_ptr);
  unsigned char GetTetraEdgeOrient(RGMTetra_ptr);
  void PrepareInputImage();

  // Utility functions
  float DistanceAtPoint(double* coords);
  float LabelAtPoint(double* coords);
  int GetIndOfMin(double, double, double, double);
  int GetIndOfMax(double, double, double, double);
  float DistanceBwPoints(double *coord0, double* coord1);
  bool SubdivInclusionTest_Resolution(RGMTetra_ptr, int);
  bool TotallyInsideBackground(RGMTetra_ptr, int);

  RGMVertex_ptr InsertVertex(double, double, double);
  RGMEdge_ptr InsertEdge(RGMVertex_ptr, RGMVertex_ptr);
  RGMTetra_ptr InsertTetraOriented(RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr);
  RGMTetra_ptr InsertTetra(unsigned char, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr);
  RGMTetra_ptr InsertTetraAllocated(int, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMEdge_ptr, RGMTetra_ptr);
  
  void RemoveBCCTetra(RGMTetra_ptr);
  void RemoveTetraAllocated(RGMTetra_ptr);
  void RemoveTetra(RGMTetra_ptr);
  void AssignLabelForTetra_Resolution_M(RGMTetra_ptr, int);
  void GetLabelToCountMapForTetra_Resolution(RGMTetra_ptr, int, std::map<int, int>&);
  void GetLabelToCountMapForTetra_Resolution2(RGMTetra_ptr, std::map<int, int>&);
  void FinalizeRedTetra(RGMTetra_ptr);
  void FinalizeGreenTetra(RGMTetra_ptr);
  void SplitEdge(RGMEdge_ptr);
  RGMEdge_ptr GetPutTmpEdge(RGMVertex_ptr, RGMVertex_ptr);
  void FindClosureRegions(bool);//true means the removal of the tetras with label 0
  void FindClosureRegionWithSpecificLabel(int);
  bool FindSubdividedRegions();
  void BuildTetraConnectivity();//build connectivity into each tetra
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryMaskTo3DAdaptiveMeshFilter.cxx"
#endif

#endif
