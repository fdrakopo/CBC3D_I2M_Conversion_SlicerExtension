/*=========================================================================
  Copyright (c) 2005 Andriy Fedorov, 
  College of William and Mary, VA and 
  Surgical Planning Lab, Harvard Medical School
  
  Adaptive 3D Tetrahedral mesh generation filter from binary mask. Algorithm 
  is described in the PhD thesis of Neil Molino (Stanford University, 2004). 
  Numerical part has been implemented by Aloys d'Aiasche. 
=========================================================================*/

#ifndef _itkBinaryMaskTo3DAdaptiveMeshFilter_txx
#define _itkBinaryMaskTo3DAdaptiveMeshFilter_txx

#include "itkBinaryMaskTo3DAdaptiveMeshFilter.h"
#include "itkNumericTraits.h"
 

namespace itk
{

template<class TInputImage, class TOutputMesh>
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::BinaryMaskTo3DAdaptiveMeshFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);
  m_SubdivInclusion = true;
  m_KeepOutside = false;
  m_NumberOfPoints = 0;
  m_NumberOfTets = 0;
  m_BCCSpacing = 10;
  m_ResampleResolution = 0;
  m_MaxAdaptiveResolution = 4;
  m_MaxResolution = 4;
  m_Fidelity = 0.8;
  m_TmpDirectory = "";
  m_InputImagePrefix = "";

  m_ThirdFaceEdgeLT[0] = -1;
  m_ThirdFaceEdgeLT[1] = -1;
  m_ThirdFaceEdgeLT[2] = 5;
  m_ThirdFaceEdgeLT[3] = 4;
  m_ThirdFaceEdgeLT[4] = 3;
  m_ThirdFaceEdgeLT[5] = 2;
  m_ThirdFaceEdgeLT[6] = -1;
  m_ThirdFaceEdgeLT[7] = -1;
  m_ThirdFaceEdgeLT[8] = 4;
  m_ThirdFaceEdgeLT[9] = 5;
  m_ThirdFaceEdgeLT[10] = 2;
  m_ThirdFaceEdgeLT[11] = 3;
  m_ThirdFaceEdgeLT[12] = 5;
  m_ThirdFaceEdgeLT[13] = 4;
  m_ThirdFaceEdgeLT[14] = -1;
  m_ThirdFaceEdgeLT[15] = -1;
  m_ThirdFaceEdgeLT[16] = 1;
  m_ThirdFaceEdgeLT[17] = 0;
  m_ThirdFaceEdgeLT[18] = 4;
  m_ThirdFaceEdgeLT[19] = 5;
  m_ThirdFaceEdgeLT[20] = -1;
  m_ThirdFaceEdgeLT[21] = -1;
  m_ThirdFaceEdgeLT[22] = 0;
  m_ThirdFaceEdgeLT[23] = 1;
  m_ThirdFaceEdgeLT[24] = 3;
  m_ThirdFaceEdgeLT[25] = 2;
  m_ThirdFaceEdgeLT[26] = 1;
  m_ThirdFaceEdgeLT[27] = 0;
  m_ThirdFaceEdgeLT[28] = -1;
  m_ThirdFaceEdgeLT[29] = -1;
  m_ThirdFaceEdgeLT[30] = 2;
  m_ThirdFaceEdgeLT[31] = 3;
  m_ThirdFaceEdgeLT[32] = 0;
  m_ThirdFaceEdgeLT[33] = 1;
  m_ThirdFaceEdgeLT[34] = -1;
  m_ThirdFaceEdgeLT[35] = -1;

  m_FaceToEdgesLT[0] = 0;
  m_FaceToEdgesLT[1] = 2;
  m_FaceToEdgesLT[2] = 5;
  m_FaceToEdgesLT[3] = 2;
  m_FaceToEdgesLT[4] = 4;
  m_FaceToEdgesLT[5] = 1;
  m_FaceToEdgesLT[6] = 4;
  m_FaceToEdgesLT[7] = 0;
  m_FaceToEdgesLT[8] = 3;
  m_FaceToEdgesLT[9] = 3;
  m_FaceToEdgesLT[10] = 5;
  m_FaceToEdgesLT[11] = 1;

  m_IncidentEdgesLT[0] = 2;
  m_IncidentEdgesLT[1] = 4;
  m_IncidentEdgesLT[2] = 2;
  m_IncidentEdgesLT[3] = 5;
  m_IncidentEdgesLT[4] = 4;
  m_IncidentEdgesLT[5] = 0;
  m_IncidentEdgesLT[6] = 4;
  m_IncidentEdgesLT[7] = 1;
  m_IncidentEdgesLT[8] = 0;
  m_IncidentEdgesLT[9] = 2;
  m_IncidentEdgesLT[10] = 0;
  m_IncidentEdgesLT[11] = 3;
  m_IncidentEdgesLT[12] = 3;
  m_IncidentEdgesLT[13] = 5;
  m_IncidentEdgesLT[14] = 3;
  m_IncidentEdgesLT[15] = 4;
  m_IncidentEdgesLT[16] = 5;
  m_IncidentEdgesLT[17] = 1;
  m_IncidentEdgesLT[18] = 5;
  m_IncidentEdgesLT[19] = 0;
  m_IncidentEdgesLT[20] = 1;
  m_IncidentEdgesLT[21] = 3;
  m_IncidentEdgesLT[22] = 1;
  m_IncidentEdgesLT[23] = 2;
}

template<class TInputImage, class TOutputMesh>
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::~BinaryMaskTo3DAdaptiveMeshFilter()
{
}

/** Generate the data */
template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::GenerateData()
 {
	if(!this->GetInput())
		itkExceptionMacro(<<"Input image missing.");

	this->Initialize();
	this->CreateMesh();

	m_CurrentResol = 0;

	// Create the helper image
	this->m_HelperImage = HelperImageType::New();
	HelperImageType::RegionType region;
	region.SetIndex(this->m_ReadyInputImage->GetLargestPossibleRegion().GetIndex());
	region.SetSize(this->m_ReadyInputImage->GetLargestPossibleRegion().GetSize());
	this->m_HelperImage->SetRegions(region);
	this->m_HelperImage->SetOrigin(this->m_ReadyInputImage->GetOrigin());
	this->m_HelperImage->SetSpacing(this->m_ReadyInputImage->GetSpacing());
	this->m_HelperImage->SetDirection(this->m_ReadyInputImage->GetDirection());
	this->m_HelperImage->Allocate();

	// Subdivide if there are still some regions needed to refined or the maximum resolution reaches.
	bool bStop = false;
	while(!bStop && m_CurrentResol < this->m_MaxResolution)
	{
		std::cout << "\nIteration " << m_CurrentResol+1 << std::endl;
		//std::cout << "Num tets : " << m_Tetras.size() << std::endl;

		// Copy all tetrahedra to a temporary list
		std::list<RGMTetra_ptr> prev_level_tetras;
		std::insert_iterator<std::list<RGMTetra_ptr > > pltI(prev_level_tetras, prev_level_tetras.begin());
		std::copy(m_Tetras.begin(), m_Tetras.end(), pltI);
		m_Tetras.clear();

		// Subdivide RED tetras (based on user criteria)
		bool new_edges_split = false;
		unsigned red_tetras_cnt = 0;
		while(prev_level_tetras.size())
		{
			RGMTetra_ptr curT;

			curT = prev_level_tetras.front();
			prev_level_tetras.pop_front();

			if(TetraUnderrefined(curT))
			{
				if(curT->parent || curT->level!=m_CurrentResol-1)
				{
					m_PendingTetras.push_back(curT);
					continue;
				}
				new_edges_split = true;
				red_tetras_cnt++;
				FinalizeRedTetra(curT);
			}
			else
				m_PendingTetras.push_back(curT);
		} // while(prev_level_size)

		//std::cout << red_tetras_cnt << " tets were Red-subdivided, "  << std::endl;

		// Enforce conformancy of the tetras neighbouring to the Red-subdivided ones. Subdivisions may incur new ones, that's why
		// there are two loops. TODO(?): bound the complexity of Red-Green subdivision procedure.
		while(new_edges_split)
		{
			new_edges_split = false;
			//std::cout << m_PendingTetras.size() << " tets are pending" << std::endl;
			assert(prev_level_tetras.empty());
			std::insert_iterator<std::list<RGMTetra_ptr > > pltI(prev_level_tetras, prev_level_tetras.begin());
			std::copy(m_PendingTetras.begin(), m_PendingTetras.end(), pltI);
			m_PendingTetras.clear();

			while(prev_level_tetras.size())
			{
				RGMTetra_ptr curT;
				unsigned split_edges_ids[6], split_edges_total;

				curT = prev_level_tetras.front();
				prev_level_tetras.pop_front();

				split_edges_total = GetTetraEdgeConf(curT);

				if(curT->parent)
				{
					// this tetra is Green: cannot subdivide, should split the parent to
					// make it Red
					// 1: disconnect all the siblings
					RGMTetra_ptr curT_parent;
					unsigned parent_edges_split, parent_child_cnt;
					unsigned child_edges_split;

					curT_parent = curT->parent;
					parent_edges_split = 0;
					parent_child_cnt = 0;

					for(unsigned i=0;i<6;i++)
						if(curT_parent->edges[i]->midv)
							parent_edges_split++;

					parent_child_cnt = curT_parent->num_greens;

					child_edges_split = 0;
					for(unsigned i=0;i<parent_child_cnt;i++)
					{
						if(&(curT_parent->greens[i])!=curT)
						{
							typename std::list<RGMTetra_ptr>::iterator gtlI = std::find(prev_level_tetras.begin(),prev_level_tetras.end(), &(curT_parent->greens[i]));
							if(gtlI==prev_level_tetras.end())
								assert(0);
							prev_level_tetras.erase(gtlI);
						}
						child_edges_split += GetTetraEdgeConf(&(curT_parent->greens[i]));
					}

					if(child_edges_split)
					{
						if(curT_parent->num_greens == 4)
						{
							for(unsigned i=1;i<=5;i+=2){
								unsigned edge_id = i;
								RGMVertex_ptr v0, v1;

								v0 = (curT_parent->greens[3]).edges[edge_id]->nodes[0];
								v1 = (curT_parent->greens[3]).edges[edge_id]->nodes[1];
								m_TmpEdgeMap[std::pair<RGMVertex_ptr,RGMVertex_ptr>(v0,v1)] = (curT_parent->greens[3]).edges[i];
							}
						} // if(curT_parent->num_greens == 4)

						for(unsigned i=0;i<parent_child_cnt;i++)
							RemoveTetraAllocated(&(curT_parent->greens[i]));
						delete [] curT_parent->greens;
						curT_parent->greens = NULL;
						FinalizeRedTetra(curT_parent);
						new_edges_split = true;
					} // if(child_edges_split)
					else
					{
						for(unsigned i=0;i<parent_child_cnt;i++)
							m_PendingTetras.push_back(&(curT_parent->greens[i]));
					} // if(child_edges_split) else
					continue;
				} // if(curT->parent)

				split_edges_total = 0;
				for(unsigned i=0;i<6;i++)
					if(curT->edges[i]->midv)
						split_edges_ids[split_edges_total++] = i;

				switch(split_edges_total)
				{
				case 0:
				case 1:
					m_PendingTetras.push_back(curT);
					break;
				case 2:
					if(split_edges_ids[0] != (split_edges_ids[1]^1))
					{
						SplitEdge(curT->edges[m_ThirdFaceEdgeLT[split_edges_ids[0]*6+split_edges_ids[1]]]);
						new_edges_split = true;
					}
					m_PendingTetras.push_back(curT);
					break;
				case 3:
					if(split_edges_ids[2] != (unsigned)m_ThirdFaceEdgeLT[split_edges_ids[0]*6+split_edges_ids[1]])
					{
						FinalizeRedTetra(curT);
						new_edges_split = true;
					}
					else
					{
						m_PendingTetras.push_back(curT);
					}
					break;
				default:
					FinalizeRedTetra(curT);
					new_edges_split = true;
				}
			} // while(prev_level_size)
		} // while(new_edges_split)

		//std::cout << m_PendingTetras.size() << " tets are pending" << std::endl;

		// Here we know that no new split edges will be introduced, and we can
		// finalize the subdivisions of Green tetrahedra.
		while(m_PendingTetras.size())
		{
			RGMTetra_ptr curT;

			curT = m_PendingTetras.front();
			m_PendingTetras.pop_front();
			if(!curT->edges_split)
				m_Tetras.push_back(curT);
			else
				FinalizeGreenTetra(curT);
		} // while(pending_tetra_cnt)

		// Initialize the helper image here
		this->m_HelperImage->FillBuffer(false);
		this->m_HelperImage->Update();

		// Find the regions that need to subdivided
		bStop = FindSubdividedRegions();
		m_CurrentResol++;
	}

	// Discard the tetrahedra, totally outside the boundary
	// Note there are two places to discard tetras. One is here, the other is in the method CreateMesh.
	//std::cout<<"\nDiscarding tetras..."<<std::endl;
	if(!m_KeepOutside)
	{
		unsigned total_tets = m_Tetras.size();
		for(unsigned i=0;i<total_tets;i++)
		{
			RGMTetra_ptr curT = m_Tetras.front();
			m_Tetras.pop_front();

			if(curT->label==0)
			{
				if(!curT->parent)
					RemoveTetra(curT);
				else
					RemoveTetraAllocated(curT);
			}
			else
				m_Tetras.push_back(curT);
		}
	}

	// Prepare the output, make sure that all terahedra are oriented
	// consistently
	std::map<RGMVertex_ptr,unsigned long> vertex2id;
	typename OutputMeshType::Pointer output_mesh = this->GetOutput();
	for(typename std::list<RGMTetra_ptr>::iterator tI=m_Tetras.begin(); tI!=m_Tetras.end(); tI++)
	{
		RGMVertex_ptr thisT_nodes[4];
		unsigned long thisT_point_ids[4];

		/* orient3d() returns a negative value if the first three points appear in
		 * counterclockwise order when viewed from the 4th point */
		if(orient3d(
				(*tI)->edges[0]->nodes[0]->coords,
				(*tI)->edges[0]->nodes[1]->coords,
				(*tI)->edges[1]->nodes[0]->coords,
				(*tI)->edges[1]->nodes[1]->coords)<0){
			thisT_nodes[0] = (*tI)->edges[0]->nodes[0];
			thisT_nodes[1] = (*tI)->edges[0]->nodes[1];
		}
		else
		{
			thisT_nodes[1] = (*tI)->edges[0]->nodes[0];
			thisT_nodes[0] = (*tI)->edges[0]->nodes[1];
		}

		thisT_nodes[2] = (*tI)->edges[1]->nodes[0];
		thisT_nodes[3] = (*tI)->edges[1]->nodes[1];

		for(unsigned i=0;i<4;i++)
		{
			if(vertex2id.find(thisT_nodes[i])==vertex2id.end())
			{
				OPointType new_point;
				new_point[0] = thisT_nodes[i]->coords[0] + m_OriginalInputOrigin[0];//restore to the original coordinate
				new_point[1] = thisT_nodes[i]->coords[1] + m_OriginalInputOrigin[1];
				new_point[2] = thisT_nodes[i]->coords[2] + m_OriginalInputOrigin[2];

				output_mesh->SetPoint(m_NumberOfPoints, new_point);
				vertex2id[thisT_nodes[i]] = m_NumberOfPoints;
				m_NumberOfPoints++;
			}
			thisT_point_ids[i] = vertex2id[thisT_nodes[i]];
		}

		TetCellAutoPointer newTet;
		newTet.TakeOwnership(new TetCell);
		newTet->SetPointIds(thisT_point_ids);
		output_mesh->SetCell(m_NumberOfTets, newTet);

		float val = (*tI)->label;
		output_mesh->SetCellData(m_NumberOfTets,val);

		m_NumberOfTets++;
	}

	std::cout << "\nNumber of Tetrahedra :" << this->m_NumberOfTets<< std::endl;
	std::cout << "Number of Vertices :" << this->m_NumberOfPoints << std::endl;

	this->GetOutput()->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());

	// Memory deallocation
	// Tets
	typename std::list<RGMTetra_ptr>::iterator tetrIter;
	for(tetrIter = this->m_Tetras.begin(); tetrIter != this->m_Tetras.end(); ++tetrIter)
	{
		if(!(*tetrIter)->parent)
			RemoveTetra(*tetrIter);
		else
			RemoveTetraAllocated(*tetrIter);
	}
	this->m_Tetras.clear();

	for(tetrIter = this->m_PendingTetras.begin(); tetrIter != this->m_PendingTetras.end(); ++tetrIter)
	{
		if(!(*tetrIter)->parent)
			RemoveTetra(*tetrIter);
		else
			RemoveTetraAllocated(*tetrIter);
	}
	this->m_PendingTetras.clear();

	// Vertices
	typename std::list<RGMVertex_ptr>::iterator vertexIter;
	for(vertexIter = this->m_Vertices.begin(); vertexIter != this->m_Vertices.end(); ++vertexIter)
		delete *vertexIter;
	this->m_Vertices.clear();

	this->m_TmpEdgeMap.clear();
	this->m_SubdivisionCriteria.clear();
 }

template<class TInputImage, class TOutputMesh>
bool
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::FindSubdividedRegions()
 {
	//std::cout << "Find the subdivided regions..." << std::endl;

	int label = *m_labelSet.begin();
	typename std::list<RGMTetra_ptr>::iterator iter;
	for (iter = m_Tetras.begin(); iter != m_Tetras.end(); ++iter)
    {
		RGMTetra_ptr curT = *iter;
		this->AssignLabelForTetra_Resolution_M(curT, m_ResampleResolution);
    }

	// Do not remove the tetras with label 0 becasue they might be included in the next iteration.
	this->FindClosureRegions(false);

	float similarityBetweenSubmeshAndMesh = 0;
	float similarityBetweenSubmeshAndRegion = 0;

	int totalNumOfVoxelsInMeshRegion = 0;
	int totalNumOfVoxelInMeshRegionWithSpecifiedLabel = 0;

	// For all tets
	for (iter = m_Tetras.begin(); iter != m_Tetras.end(); ++iter)
	{
		RGMTetra_ptr curT = *iter;
		if(curT->label != label)
			continue;

		std::map<int, int>::iterator it;
		std::map<int, int> labelToCountMap;
		this->GetLabelToCountMapForTetra_Resolution2(curT, labelToCountMap);

		for(it=labelToCountMap.begin(); it!=labelToCountMap.end(); it++)
		{
			totalNumOfVoxelsInMeshRegion += (*it).second;
			if((*it).first == label)
				totalNumOfVoxelInMeshRegionWithSpecifiedLabel += (*it).second;
		}
	}

	// three numbers:
	// 1. the total number of the voxels in the mesh region
	// 2. the total number of the voxels in the mesh region with current region label
	// 3. the total number of the voxels in the real region (label image)
	if (totalNumOfVoxelsInMeshRegion == 0)
		similarityBetweenSubmeshAndMesh = 0;
	else
		similarityBetweenSubmeshAndMesh = (float)totalNumOfVoxelInMeshRegionWithSpecifiedLabel/totalNumOfVoxelsInMeshRegion;

	similarityBetweenSubmeshAndRegion = (float)totalNumOfVoxelInMeshRegionWithSpecifiedLabel/m_labelToCount[label];

	//std::cout << "\nTissue with label: " << label << std::endl;
	//std::cout << "Number of common voxels between tissue and mesh (S1∩S2) : " << totalNumOfVoxelInMeshRegionWithSpecifiedLabel<<std::endl;
	//std::cout << "Number of voxels in the mesh region (S1) : " << totalNumOfVoxelsInMeshRegion << std::endl;
	//std::cout << "Number of voxels in the tissue region (S2) : " << m_labelToCount[label]<< std::endl;

	// checks
	if(similarityBetweenSubmeshAndMesh > 1.)
		itkExceptionMacro("Similarity Between Submesh And Mesh is " << similarityBetweenSubmeshAndMesh << ". Should be less than 1.");
	if(similarityBetweenSubmeshAndRegion > 1.)
		itkExceptionMacro("Similarity Between Submesh And Region is " << similarityBetweenSubmeshAndRegion << ". Should be less than 1.");

	if(similarityBetweenSubmeshAndMesh <= m_Fidelity)
		std::cout << "Similarity between common region and mesh region: " << similarityBetweenSubmeshAndMesh << " < " << m_Fidelity << std::endl;
	else
		std::cout << "Similarity between common region and mesh region: " << similarityBetweenSubmeshAndMesh << " > " << m_Fidelity << std::endl;

	if(similarityBetweenSubmeshAndRegion <= m_Fidelity)
		std::cout << "Similarity between common region and tissue region: " << similarityBetweenSubmeshAndRegion << " < " << m_Fidelity << std::endl;
	else
		std::cout << "Similarity between common region and tissue region: " << similarityBetweenSubmeshAndRegion << " > " << m_Fidelity << std::endl;

	if(similarityBetweenSubmeshAndMesh > m_Fidelity && similarityBetweenSubmeshAndRegion > m_Fidelity)
	{
		//std::cout<<"Tissue region with label " << label << " is finalized (resolution = " << m_CurrentResol << ")" << std::endl;
		return true;
	}
	else
	{
		//std::cout<<"Tissue region with label " << label << " is NOT finalized (resolution = " << m_CurrentResol << ")" << std::endl;
		return false;
	}
}

//find surfaces and vertices. The surfaces will define multiply regions, each of which is characterized by
//for each tetra, at most one boundary face. If a tetra has more than one boundary faces, this tetra needs
//to be further divided (several different boundary faces) or belong to another region (several save boundary faces).
template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::FindClosureRegions(bool isBackgroundTetraDiscarded)
 {
	//1. save connectivity into each tetra
	this->BuildTetraConnectivity();

	//2. the label loop is used to find the closure surface. The labels will be redistributed in the mesh
	this->m_finalizedLabelSet.clear();
	for (std::set<int>::iterator it = m_labelSet.begin(); it != m_labelSet.end(); it++)
	{
		if((*it) != 0)
		{
			this->FindClosureRegionWithSpecificLabel(*it);
			this->m_finalizedLabelSet.insert((*it));
		}
	}

	//3. Discard the tetrahedra with label 0
	if( !isBackgroundTetraDiscarded)
		return;

	unsigned int total_tets = this->m_Tetras.size();
	RGMTetra_ptr curTet;
	int i;
	for(i=0; i<(int)total_tets; i++)
	{
		curTet = this->m_Tetras.front();
		m_Tetras.pop_front();

		if(curTet->label==0 )
		{
			if(!curTet->parent)
				RemoveTetra(curTet);
			else
				RemoveTetraAllocated(curTet);
		}
		else
			this->m_Tetras.push_back(curTet);
	}
}

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::BuildTetraConnectivity()
 {
	// Walk through the input mesh, extract surface vertices

	typedef typename std::map< TetraRGFace,FaceInfo* > Face2InfoMap;
	typedef typename std::map< TetraRGFace,FaceInfo* >::iterator Face2InfoMapIterator;

	Face2InfoMap face2info;
	Face2InfoMapIterator face2infoIterator;
	std::set<int> labelSet;
	unsigned int j;

	const unsigned int tet_face_LUT[16] =
			{ 0, 1, 2, 3,
			0, 1, 3, 2,
			0, 2, 3, 1,
			1, 2, 3, 0 };

	//std::cout << "Extracting surface vertices..." << std::endl;

	RGMVertex_ptr vertices[4];
	RGMTetra_ptr curTet;
	int label;

	//1. put tetra into different list and (for the simplicity, remove this functionality)
	//build share information for each face and build label information
	typename std::list<RGMTetra_ptr>::iterator iter;
	for (iter = this->m_Tetras.begin(); iter != this->m_Tetras.end(); ++iter)
	{
		curTet = *iter;

		//collect label information
		label = curTet->label;
		labelSet.insert(label);

		//build the "share information" for each face
		vertices[0] = curTet->edges[0]->nodes[0];
		vertices[1] = curTet->edges[0]->nodes[1];
		vertices[2] = curTet->edges[1]->nodes[0];
		vertices[3] = curTet->edges[1]->nodes[1];

		for(j=0;j<4;j++)
		{
			TetraRGFace thisFace(vertices[tet_face_LUT[j*4]],
					vertices[tet_face_LUT[j*4+1]],
					vertices[tet_face_LUT[j*4+2]],
					vertices[tet_face_LUT[j*4+3]]);

			face2infoIterator = face2info.find(thisFace);
			if(face2infoIterator!=face2info.end())
			{
				FaceInfo *fi = (*face2infoIterator).second;
				//fi->n++;//n should be 2
				fi->ptr1 = curTet;
			}
			else
			{
				FaceInfo *fi = new FaceInfo;
				//fi->n=1;
				fi->ptr0 = curTet;
				fi->ptr1 = NULL;
				face2info[thisFace] = fi;
			}
		}
	}

	//2. based on the share information of each face, build the share information for each tetra
	for (iter = this->m_Tetras.begin(); iter != this->m_Tetras.end(); ++iter)
	{
		curTet = *iter;
		vertices[0] = curTet->edges[0]->nodes[0];
		vertices[1] = curTet->edges[0]->nodes[1];
		vertices[2] = curTet->edges[1]->nodes[0];
		vertices[3] = curTet->edges[1]->nodes[1];

		//curTet->num_sharedTetra = 0;
		for(j=0;j<4;j++)
		{
			TetraRGFace thisFace(vertices[tet_face_LUT[j*4]],
					vertices[tet_face_LUT[j*4+1]],
					vertices[tet_face_LUT[j*4+2]],
					vertices[tet_face_LUT[j*4+3]]);

			face2infoIterator = face2info.find(thisFace);

			if(face2infoIterator!=face2info.end())
			{
				FaceInfo *fi = (*face2infoIterator).second;
				if(fi->ptr1==NULL)//prt0 must be non null
					curTet->s[j] = NULL;
				else//n==2
				{
					if(curTet==fi->ptr0)
						curTet->s[j]=fi->ptr1;
					else
						curTet->s[j]=fi->ptr0;
					//curTet->num_sharedTetra++;
				}
			}
			else
			{
				//std::cout << "Should find this face!!!" << std::endl;
			}
		}
	}
	//release FaceInfo in face2info
	for(face2infoIterator=face2info.begin(); face2infoIterator!=face2info.end(); face2infoIterator++)
	{
		FaceInfo *fi = (*face2infoIterator).second;
		delete fi;
	}
}

//The principle to find the closure region
//1. no connectivity assumption for each tetra
//2. put boudary tetras into to neighboring tissues or background as can as possible???
template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::FindClosureRegionWithSpecificLabel(int _label)
 {
	typedef typename std::map<int, int> Label2CountMap;
	typedef typename std::map<int, int>::iterator Label2CountMapIterator;

	int i, j, neighboringLabel;
	bool isChange = true;

	RGMTetra_ptr curTet;

	unsigned total_tets = m_Tetras.size();

	while(isChange)
	{
		isChange = false;

		for(i=0;i<(int)total_tets;i++)
		{

			Label2CountMap label2Count;
			Label2CountMapIterator label2CountIterator;


			curTet = m_Tetras.front();

			m_Tetras.pop_front();

			if(curTet->label != _label)
			{
				m_Tetras.push_back(curTet);
				continue;
			}

			//obtain label to count for the tetra
			for(j=0; j<4; j++)
			{
				if(curTet->s[j]==NULL)
					neighboringLabel=0;
				else
					neighboringLabel = curTet->s[j]->label;
				label2CountIterator = label2Count.find(neighboringLabel);
				if(label2CountIterator!=label2Count.end())
					label2Count[neighboringLabel]++;//good way to change the value in the map?
				else
					label2Count[neighboringLabel]=1;
			}

			//the sum of the count must be 4
			//the number of the shared region cannot exceed 2, otherwise, this tetra needs to be subdivided.
			//The method performs with a greedy fashion. For each region, put the tetras, which do not satisfy the requirement
			//for the closure region into other regions.
			switch(label2Count.size())
			{
			case 0:{
				std::cout << "The number of the shared tetra with the samel label cannot be zero!" << std::endl;
				break;
			}
			case 1://four shared tetras have the same label
			{
				label2CountIterator=label2Count.begin();
				neighboringLabel=(*label2CountIterator).first;
				if(curTet->label != neighboringLabel)
				{
					//std::cout<<"The tetra with label " << curTet->label <<" is a hole with neighboring: "<< neighboringLabel<<std::endl;
					curTet->label=neighboringLabel;
					isChange =true;
				}
				break;
			}
			case 2://there are two cases: the count can be (1, 3) or (2, 2)
			{

				label2CountIterator=label2Count.begin();
				int neighboringLabel1 = (*label2CountIterator).first;
				int count1 = (*label2CountIterator).second;
				label2CountIterator++;
				int neighboringLabel2 = (*label2CountIterator).first;
				int count2 = (*label2CountIterator).second;

				if(count1==count2)
				{
					if(curTet->label!=neighboringLabel1 && curTet->label!= neighboringLabel2)
					{
						if(m_finalizedLabelSet.find(neighboringLabel1) == m_finalizedLabelSet.end())
							curTet->label = neighboringLabel1;
						else
							if(m_finalizedLabelSet.find(neighboringLabel2) == m_finalizedLabelSet.end())
								curTet->label = neighboringLabel2;
							else
								curTet->label = neighboringLabel1;//for simplicity, do not obey to principle 2 here

						isChange =true;

					}
					else
					{
						if(neighboringLabel1 != curTet->label)
							if(m_finalizedLabelSet.find(neighboringLabel1) == m_finalizedLabelSet.end())
							{
								curTet->label = neighboringLabel1;
								isChange = true;
							}
							else{}
						else
							if(m_finalizedLabelSet.find(neighboringLabel2) == m_finalizedLabelSet.end())
							{
								curTet->label = neighboringLabel2;
								isChange = true;
							}
							else{}
					}
				}
				else// (1,3)//in this case, do not consider whether the neighboring regions are finalized.
				{
					if(count1==3)
						if(curTet->label != neighboringLabel1)//3 exposed faces, must remove, so do not consider the finalized labels
						{
							curTet->label = neighboringLabel1;
							isChange = true;
						}
						else//three faces are in the connected region and one is outside
						{}
					else//count2==3
						if(curTet->label != neighboringLabel2)//3 exposed faces
						{
							curTet->label = neighboringLabel2;
							isChange = true;
						}
						else//three faces are in the connected region and one is outside
						{}

				}

				break;

			}
			case 3:// (1, 1, 2). Each tetra only can be connected with at most one other region
			{
				//std::cout << "This tetra needs to be subdivided!" << std::endl;
				//put it into the non-finalized region if available.
				label2CountIterator=label2Count.begin();
				int neighboringLabel1 = (*label2CountIterator).first;
				label2CountIterator++;
				int neighboringLabel2 = (*label2CountIterator).first;
				label2CountIterator++;
				int neighboringLabel3 = (*label2CountIterator).first;


				if(curTet->label!=neighboringLabel1 && curTet->label!=neighboringLabel2 && curTet->label!=neighboringLabel3)
				{
					if(m_finalizedLabelSet.find(neighboringLabel1) == m_finalizedLabelSet.end())
						curTet->label = neighboringLabel1;
					else
						if(m_finalizedLabelSet.find(neighboringLabel2) == m_finalizedLabelSet.end())
							curTet->label = neighboringLabel2;
						else
							if(m_finalizedLabelSet.find(neighboringLabel3) == m_finalizedLabelSet.end())
								curTet->label = neighboringLabel3;
							else
								curTet->label = neighboringLabel1;//for simplicity, do not obey to principle 2 here

					isChange =true;

				}
				else
				{
					for(label2CountIterator=label2Count.begin(); label2CountIterator!=label2Count.end(); label2CountIterator++)
					{
						neighboringLabel=(*label2CountIterator).first;
						if(curTet->label == neighboringLabel)
							continue;

						if(m_finalizedLabelSet.find(neighboringLabel) == m_finalizedLabelSet.end())
						{
							curTet->label = neighboringLabel;
							isChange = true;

							break;
						}
					}
				}
				break;
			}
			case 4:
			{
				//put it into the non-finalized region if available.
				label2CountIterator=label2Count.begin();
				int neighboringLabel1 = (*label2CountIterator).first;
				label2CountIterator++;
				int neighboringLabel2 = (*label2CountIterator).first;
				label2CountIterator++;
				int neighboringLabel3 = (*label2CountIterator).first;
				label2CountIterator++;
				int neighboringLabel4 = (*label2CountIterator).first;

				if(curTet->label!=neighboringLabel1 && curTet->label!=neighboringLabel2 && curTet->label!=neighboringLabel3 && curTet->label!=neighboringLabel4)
				{
					if(m_finalizedLabelSet.find(neighboringLabel1) == m_finalizedLabelSet.end())
						curTet->label = neighboringLabel1;
					else
						if(m_finalizedLabelSet.find(neighboringLabel2) == m_finalizedLabelSet.end())
							curTet->label = neighboringLabel2;
						else
							if(m_finalizedLabelSet.find(neighboringLabel3) == m_finalizedLabelSet.end())
								curTet->label = neighboringLabel3;
							else
								if(m_finalizedLabelSet.find(neighboringLabel4) == m_finalizedLabelSet.end())
									curTet->label = neighboringLabel4;
								else
									curTet->label = neighboringLabel1;//for simplicity, do not obey to principle 2 here

					isChange =true;
				}
				else
				{
					for(label2CountIterator=label2Count.begin(); label2CountIterator!=label2Count.end(); label2CountIterator++)
					{
						neighboringLabel=(*label2CountIterator).first;
						if(curTet->label == neighboringLabel)
							continue;

						if(m_finalizedLabelSet.find(neighboringLabel) == m_finalizedLabelSet.end())
						{
							curTet->label = neighboringLabel;
							isChange = true;

							break;

						}
					}
				}
				//std::cout << "This tetra needs to be subdivided!" << std::endl;
				break;
			}

			default:
			{}
			}
			m_Tetras.push_back(curTet);
		}
	}
}

/** Initialize the class fields */
template<class TInputImage, class TOutputMesh>
bool
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::Initialize()
 {
	//std::cout << "Initializing..." << std::endl;

	// Initialize non-constant class members
	m_Interpolator = InterpolatorType::New();

	// Set the input image
	this->m_InputImage = static_cast<InputImageType*>(this->ProcessObject::GetInput(0));

	// Cast the input image to the internal format
	if(m_InputImagePrefix.size())
	{
		// resampled to unit-voxel input image
		try
		{
			InternalImageReaderType::Pointer reader = InternalImageReaderType::New();
			reader->SetFileName((m_TmpDirectory+m_InputImagePrefix+"Ready.mhd").c_str());
			reader->Update();
			this->m_ReadyInputImage = reader->GetOutput();
		}
		catch(ExceptionObject &e)
		{
			InternalImageWriterType::Pointer writer = InternalImageWriterType::New();
			PrepareInputImage();
			writer->SetFileName((m_TmpDirectory+m_InputImagePrefix+"Ready.mhd").c_str());
			writer->SetInput(this->m_ReadyInputImage);
			try
			{
				writer->Update();
			}
			catch(ExceptionObject &eo)
			{
				// this should not disrupt the overall execution
			}
		}
	}
	else
	{
		PrepareInputImage();
	}

	typename InternalImageType::SizeType input_size = this->m_ReadyInputImage->GetLargestPossibleRegion().GetSize();
	this->m_dimX = input_size[0];
	this->m_dimY = input_size[1];
	this->m_dimZ = input_size[2];

	m_OriginalInputOrigin = m_InputImage->GetOrigin();//restore original origin
	m_InputOrigin[0]=0.0,  m_InputOrigin[1]=0.0,  m_InputOrigin[2]=0.0;//use the regularized image (origin=0 and spacing=1) for meshing

	// Initialize the interpolator
	m_Interpolator->SetInputImage(m_ReadyInputImage);

	//Initialize the label set
	float *imageData = (float*)m_ReadyInputImage->GetBufferPointer();
	std::map<int, int>::iterator it;
	int pixel;
	for(int i=0; i<m_dimX*m_dimY*m_dimZ; i++)
	{
		//used for the finding of the closure region
		pixel = (int)imageData[i];
		m_labelSet.insert(pixel);

		//used for the subdivision
		if(m_labelToCount.find(pixel) != m_labelToCount.end())
			m_labelToCount[pixel]++;
		else
			m_labelToCount[pixel] = 1;
	}

	int numLabels = m_labelToCount.size();
	if(numLabels != 2)
		itkExceptionMacro("Invalid input image! " << numLabels << " labels detected (including the background). "
				"Currently the method supports only 2 labels (background and tissue).");

	// Remove the background label
	m_labelSet.erase(0);

	// Find out the appropriate value of the BCC spacing
	m_BCCSpacing = FindBCCSpacing();
	//std::cout << "Initialization done" << std::endl;
	return 0;
}

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::CreateMesh()
 {
	m_CurrentResol = 0;
	if (this->GetNumberOfInputs() < 1)
		itkExceptionMacro(<< "Invalid number of inputs.");

	typename InputImageType::ConstPointer m_InputImage = static_cast<const InputImageType * >(this->ProcessObject::GetInput(0) );
	if(!this->m_InputImage)
		itkExceptionMacro(<< "Input Image is Null.");

	this->CreateBCC();
	//std::cout << "BCC created" << std::endl;

	// If requested, need to remove the elements which are guaranteed to be outside the object
	unsigned long total_tets = 0;
	if(!m_KeepOutside)
	{
		// "The tetrahedron is guaranteed not to intersect the interface
		// if the minimum valued of |phi| at a node is larger than the longest
		// edge length" [Molino]
		// Of course, a tetrahedron is not discarded if it has both negative and
		// positive values of distance at the vertices.
		total_tets = m_Tetras.size();
		//std::cout << "Number of tets before removals: " << total_tets << std::endl;
		unsigned long i;

		for(i=0; i<total_tets; i++)
		{
			RGMTetra_ptr thisT = m_Tetras.front();
			m_Tetras.pop_front();

			// Check if the tetrahedron is totally inside the background
			if(TotallyInsideBackground(thisT, m_ResampleResolution))
				RemoveBCCTetra(thisT);
			else
				m_Tetras.push_back(thisT);
		}
	}
	total_tets = m_Tetras.size();
	//std::cout << "Number of tets after removals : " << total_tets << std::endl;
}

/** PrintSelf */
template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);

  /*
  os << indent
     << "ObjectValue: " 
     << static_cast<NumericTraits<unsigned char>::PrintType>(m_ObjectValue)
     << std::endl;

  os << indent
     << "NumberOfNodes: "
     << m_NumberOfNodes
     << std::endl;

  os << indent
     << "NumberOfCells: "
     << m_NumberOfCells
     << std::endl;
     */
}

/* Returns the interpolated value of distance at the specified coordinate.
 * Will perform the inside test.
 */
template<class TInputImage, class TOutputMesh>
float
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::DistanceAtPoint(double *coords){

  float distance = 0;
  typename InterpolatorType::ContinuousIndexType input_index;

  input_index[0] = coords[0];
  input_index[1] = coords[1];
  input_index[2] = coords[2];

  if(m_Interpolator->IsInsideBuffer(input_index))
    distance = (float)m_Interpolator->EvaluateAtContinuousIndex(input_index);
  else {
    std::cerr << "DistanceAtPoint(): Point [" << coords[0] << ", " << coords[1] << ", " 
      << coords[2] << "] is outside the image boundaries" << std::endl;
    assert(0);
  }
  return distance;
}

/* Returns the interpolated value (nearest neighborhood) of label at the specified coordinate.
 * Will perform the inside test.
 */
template<class TInputImage, class TOutputMesh>
float
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::LabelAtPoint(double *coords)
 {
	float label = 0;
	typename InterpolatorType::ContinuousIndexType input_index;

	//The coords can be used  as the index, because the image has been regularized.
	input_index[0] = coords[0];
	input_index[1] = coords[1];
	input_index[2] = coords[2];

	if(m_Interpolator->IsInsideBuffer(input_index))
	{
		label = (float)m_Interpolator->EvaluateAtContinuousIndex(input_index);
	}
	else
	{
		std::cerr << "LabelAtPoint(): Point [" << coords[0] << ", " << coords[1] << ", "
				<< coords[2] << "] is outside the image boundaries" << std::endl;
		assert(0);
	}
	return label;
}

//Assign the label to count map for the tetra with specific resolution
//This method is used for the determinement of the subdivision and the assignment of the label.
template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::GetLabelToCountMapForTetra_Resolution(RGMTetra_ptr thisT, int resolution, std::map<int, int> &labelToCountMap)
 {
	int k,l,m;
	int Indices[3][2];
	float x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,a,b,c,d;
	unsigned int check = 0;
	double v[4][3];
	std::map<int, int>::iterator it;

	float label=0;
	double coords[3];

	labelToCountMap.clear();

	RGMVertex_ptr vertices[4];
	vertices[0] = thisT->edges[0]->nodes[0];
	vertices[1] = thisT->edges[0]->nodes[1];
	vertices[2] = thisT->edges[1]->nodes[0];
	vertices[3] = thisT->edges[1]->nodes[1];

	int scaling = 1 << resolution;

	//change to the coordinates under new resolution
	for(unsigned int i=0; i<4; i++)
		for(unsigned int j=0; j<3; j++)
			v[i][j]=vertices[i]->coords[j]*scaling;


	for (unsigned int i=0;i<3;i++) /* for x,y and z coordinates of tetrahedron */
	{
		Indices[i][0] = GetIndOfMin(v[0][i],v[1][i],v[2][i],v[3][i]);
		Indices[i][1] = GetIndOfMax(v[0][i],v[1][i],v[2][i],v[3][i]);
	}

	assert( (Indices[0][0] >=0 ) && (Indices[0][1] < this->m_dimX * scaling) );
	assert( (Indices[1][0] >=0 ) && (Indices[1][1] < this->m_dimY * scaling) );
	assert( (Indices[2][0] >=0 ) && (Indices[2][1] < this->m_dimZ * scaling) );

	////	std::cout << v[0][0] << " " << v[0][1] << " " << v[0][2] << std::endl;
	////	std::cout << v[1][0] << " " << v[1][1] << " " << v[1][2] << std::endl;
	////	std::cout << v[2][0] << " " << v[2][1] << " " << v[2][2] << std::endl;
	////	std::cout << v[3][0] << " " << v[3][1] << " " << v[3][2] << std::endl;
	////	std::cout << " " << std::endl;
	//	int numPoints = 0;

	/* for every point k,l,m inside the tet's cube*/
	for (k=Indices[0][0];k<(Indices[0][1]+1);k++)
	{
		for (l=Indices[1][0];l<(Indices[1][1]+1);l++)
		{
			for (m=Indices[2][0];m<(Indices[2][1]+1);m++)
			{
				//				std::cout << k << " " << l << " " << m << std::endl;
				//				numPoints++;
				check = 0;
				for (unsigned int i=0; i<4; i++)
				{
					/* x1,x2,x3 = plane, x4 = 4th point */
					x1 = v[i%4][0];
					y1 = v[i%4][1];
					z1 = v[i%4][2];

					x2 = v[(i+1)%4][0];
					y2 = v[(i+1)%4][1];
					z2 = v[(i+1)%4][2];

					x3 = v[(i+2)%4][0];
					y3 = v[(i+2)%4][1];
					z3 = v[(i+2)%4][2];

					x4 = v[(i+3)%4][0];
					y4 = v[(i+3)%4][1];
					z4 = v[(i+3)%4][2];

					/* Coefficients of the plane's equation */
					a = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1) ;
					b = (z2-z1)*(x3-x1) - (z3-z1)*(x2-x1) ;
					c = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) ;
					d = -x1*a - y1*b -z1*c ;

					/*Check if point is on the same side as corresponding plane */
					if (((a*x4 + b*y4 + c*z4 + d)*(a*k + b*l + c*m + d))>=0)
						check++;
				}

				// The point is inside
				if (check == 4)
				{
					//get coords with the original resolution
					coords[0] = (double)k/scaling;
					coords[1] = (double)l/scaling;
					coords[2] = (double)m/scaling;

					label = this->LabelAtPoint(coords);

					if(labelToCountMap.find((int)label)!=labelToCountMap.end())
						labelToCountMap[(int)label]++;
					else
						labelToCountMap[(int)label] = 1;
				}
			}
		}
	}
	//std::cout << "Num points : " << numPoints << std::endl;
 }

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::GetLabelToCountMapForTetra_Resolution2(RGMTetra_ptr thisT,std::map<int, int> &labelToCountMap)
 {
	int k,l,m;
	int Indices[3][2];
	float x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,a,b,c,d;
	unsigned int count = 0;
	double v[4][3];
	float label=0;
	labelToCountMap.clear();

	RGMVertex_ptr vertices[4];
	vertices[0] = thisT->edges[0]->nodes[0];
	vertices[1] = thisT->edges[0]->nodes[1];
	vertices[2] = thisT->edges[1]->nodes[0];
	vertices[3] = thisT->edges[1]->nodes[1];

	for(unsigned int i=0; i<4; i++)
		for(unsigned int j=0; j<3; j++)
			v[i][j] = vertices[i]->coords[j];

	for (unsigned int i=0;i<3;i++)
	{
		Indices[i][0] = GetIndOfMin(v[0][i],v[1][i],v[2][i],v[3][i]);
		Indices[i][1] = GetIndOfMax(v[0][i],v[1][i],v[2][i],v[3][i]);
	}

	HelperImageType::IndexType index;
	HelperImageType::RegionType region = this->m_HelperImage->GetLargestPossibleRegion();

	// For every point k,l,m inside the tet's cube
	for (k = Indices[0][0]; k < (Indices[0][1]+1); k++)
	{
		// check
		if(k < 0)
			itkExceptionMacro("index k is " << k << ". Should be >= 0");

		index[0] = k;

		for (l = Indices[1][0]; l < (Indices[1][1]+1); l++)
		{
			// check
			if(l < 0)
				itkExceptionMacro("index l is " << l << ". Should be >= 0");

			index[1] = l;

			for (m = Indices[2][0]; m < (Indices[2][1]+1); m++)
			{
				// check
				if(m < 0)
					itkExceptionMacro("index m is " << m << ". Should be >= 0");

				index[2] = m;

				// Check if the index is inside the image region
				if(!region.IsInside(index))
					continue;

				// Check if the index has already set to true
				if(this->m_HelperImage->GetPixel(index))
					continue;

				count = 0;
				for (unsigned int i=0; i<4; i++)
				{
					/* x1,x2,x3 = plane, x4 = 4th point */
					x1 = v[i%4][0];
					y1 = v[i%4][1];
					z1 = v[i%4][2];

					x2 = v[(i+1)%4][0];
					y2 = v[(i+1)%4][1];
					z2 = v[(i+1)%4][2];

					x3 = v[(i+2)%4][0];
					y3 = v[(i+2)%4][1];
					z3 = v[(i+2)%4][2];

					x4 = v[(i+3)%4][0];
					y4 = v[(i+3)%4][1];
					z4 = v[(i+3)%4][2];

					/* Coefficients of the plane's equation */
					a = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1) ;
					b = (z2-z1)*(x3-x1) - (z3-z1)*(x2-x1) ;
					c = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) ;
					d = -x1*a - y1*b -z1*c ;

					/*Check if point is on the same side as corresponding plane */
					if (((a*x4 + b*y4 + c*z4 + d)*(a*k + b*l + c*m + d))>=0)
						count++;
				}

				// The point is inside
				if (count == 4)
				{
					// Get the pixel value
					// The index is guaranteed to be inside the image here because the helper image and
					// the m_ReadyInputImage have same region.
					label = this->m_ReadyInputImage->GetPixel(index);

					// Check
					if(this->m_HelperImage->GetPixel(index))
						itkExceptionMacro("The index has already set to true!");

					if(labelToCountMap.find((int)label) != labelToCountMap.end())
						labelToCountMap[(int)label]++;
					else
						labelToCountMap[(int)label] = 1;

					// The index was checked
					this->m_HelperImage->SetPixel(index,true);
				}
			}
		}
	}
 }

//Assign the label with maximum number of the voxel in the tetra to the tetra
template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::AssignLabelForTetra_Resolution_M(RGMTetra_ptr thisT, int resolution)
 {
	//find max label
	int MAXCount = -1000;
	int label=0;
	std::map<int, int> labelToCountMap;
	std::map<int, int>::iterator it;
	int adaptive_resolution = resolution;

	if(adaptive_resolution > (int)this->m_MaxAdaptiveResolution)
		itkExceptionMacro("Adaptive resolution (=" << adaptive_resolution << ") is larger than the Max resolution (=" << this->m_MaxAdaptiveResolution << ")");

	while(adaptive_resolution >= 0 && adaptive_resolution <= (int)this->m_MaxAdaptiveResolution)
	{
		this->GetLabelToCountMapForTetra_Resolution(thisT, adaptive_resolution, labelToCountMap);
		if(labelToCountMap.size() > 0)
			break;

		adaptive_resolution++;
	}

	if(labelToCountMap.size() == 0)
		itkExceptionMacro("There are no voxels inside the tetra! Need to specify larger Maximum resolution!");

	for(it=labelToCountMap.begin();it!=labelToCountMap.end();it++)
	{
		if((*it).second > MAXCount)
		{
			MAXCount= (*it).second;
			label = (*it).first;
		}
	}
	thisT->label = label;
 }

/********************************************
 *	get the min of the four value
 ********************************************/
template<class TInputImage, class TOutputMesh>
int
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::GetIndOfMin(double a,double b,double c,double d)
  {
    if ((a<=b)&&(a<=c)&&(a<=d)) {
      return (int)(a);
    }
    if ((b<=a)&&(b<=c)&&(b<=d)){
      return (int)(b);
    }
    if ((c<=d)&&(c<=a)&&(c<=b)){
      return (int)(c);
    }
    if ((d<=a)&&(d<=b)&&(d<=c)){
      return (int)(d);
    }

	return 0;//only for correctly compiling
  }//end of GetIndOfMin


/********************************************
 *	get the max of the four value
 ********************************************/
template<class TInputImage, class TOutputMesh>
int
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::GetIndOfMax(double a,double b,double c,double d)
  {
    if ((a>=b)&&(a>=c)&&(a>=d))
      return ((int)floor(a));
    if ((b>=a)&&(b>=c)&&(b>=d))
      return ((int)floor(b));
    if ((c>=d)&&(c>=a)&&(c>=b))
      return ((int)floor(c)); 
    if ((d>=a)&&(d>=b)&&(d>=c))
      return ((int)(d));

	return 0;//only for correctly compiling
  }//end of GetIndOfMax

/* Computes the distance between two points in space */
template<class TInputImage, class TOutputMesh>
float
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::DistanceBwPoints(double *c0, double *c1){
  return sqrt((c0[0]-c1[0])*(c0[0]-c1[0])+(c0[1]-c1[1])*(c0[1]-c1[1])+
    (c0[2]-c1[2])*(c0[2]-c1[2]));
}

/* Return true if the tetraheron crosses the surface and has to 
 * be subdivided. */
template<class TInputImage, class TOutputMesh>
bool
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::SubdivInclusionTest_Resolution(RGMTetra_ptr thisT, int resolution)
 {
	std::map<int, int> labelToCountMap;
	int adaptive_resolution = resolution;

	if(adaptive_resolution > (int)this->m_MaxAdaptiveResolution)
		itkExceptionMacro("Adaptive resolution (=" << adaptive_resolution << ") is larger than the Max resolution (=" << this->m_MaxAdaptiveResolution << ")");

	while(adaptive_resolution >= 0 && adaptive_resolution <= (int)this->m_MaxAdaptiveResolution)
	{
		this->GetLabelToCountMapForTetra_Resolution(thisT, adaptive_resolution, labelToCountMap);
		if(labelToCountMap.size() > 0)
			break;

		adaptive_resolution++;
	}

	// Checks
	if(labelToCountMap.size() == 0)
		itkExceptionMacro("There are no voxels inside the tetra! Need to specify larger Maximum resolution!");
	if(labelToCountMap.size() > 2)
		itkExceptionMacro("Invalid number of labels inside the tetra!");

	// The terahedron will not subdivided (one label).
	if(labelToCountMap.size() == 1)
		return false;
	// The terahedron will subdivided (two labels).
	else
		return true;
}

/* Return true if all the voxels contained in this tetra are background */
template<class TInputImage, class TOutputMesh>
bool
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::TotallyInsideBackground(RGMTetra_ptr thisT, int resolution)
 {
  std::map<int, int> labelToCountMap;
  std::map<int, int>::iterator it;

  this->GetLabelToCountMapForTetra_Resolution(thisT, resolution, labelToCountMap);

  if(labelToCountMap.size() == 0)
	  itkExceptionMacro("There is no voxel in the tetra! Need to specify a larger input Resample Resolution (1 or 2)!");

  if(labelToCountMap.size() == 1 && ((*labelToCountMap.begin()).first == 0))//this tetra only includes background voxels
	  return true;
  else
	  return false;
}

/* Generates the initial BCC lattice */
template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::CreateBCC(){
  // TODO: make the BCC conformant with the image boundaries

  // The outer BCC is 1 spacing distance larger than the inner BCC in all
  // dimensions. The number of tetrahedra in level 0 can be calculated 
  // Start traversal of the vertices of outer BCC, the first vertex has
  // coordinates (m_x0+spacing,m_y0+spacing,m_z0+spacing). 

  // Generate the BCC vertices
  double dx, dy, dz;
  unsigned i, j, k;
  unsigned dsi, dsj, dsk;
  double half_spacing = this->m_BCCSpacing/2.0;
  double m_x0, m_y0, m_z0;
  
  m_x0 = m_InputOrigin[0];
  m_y0 = m_InputOrigin[1];
  m_z0 = m_InputOrigin[2];

  dsi = (int) ((this->m_dimX-1)/this->m_BCCSpacing) + 1;
  dsj = (int) ((this->m_dimY-1)/this->m_BCCSpacing) + 1;
  dsk = (int) ((this->m_dimZ-1)/this->m_BCCSpacing) + 1;

  RGMVertex_ptr* BCC = new RGMVertex_ptr[(dsi+1)*(dsj+1)*(dsk+1)];
  RGMVertex_ptr* iBCC = new RGMVertex_ptr[dsi*dsj*dsk];

  for(dz=m_z0,k=0;k<dsk;dz+=m_BCCSpacing,k++){
    for(dy=m_y0,j=0;j<dsj;dy+=m_BCCSpacing,j++){
      for(dx=m_x0,i=0;i<dsi;dx+=m_BCCSpacing,i++){
        BCC[k*(dsi+1)*(dsj+1)+j*(dsi+1)+i] = InsertVertex(dx, dy, dz);
        if(dx!=m_x0+m_dimX && dy!=m_y0+m_dimY && dz!=m_z0+m_dimZ){
          iBCC[k*dsi*dsj+j*dsi+i] = InsertVertex(dx+half_spacing, dy+half_spacing, dz+half_spacing);
        }
      }
    }
  }

  int ii, ik, ij;
  RGMEdge_ptr *BCC_edges = new RGMEdge_ptr[dsj*dsk*dsi*3];
  RGMEdge_ptr *iBCC_edges = new RGMEdge_ptr[dsj*dsk*dsi*3];
  RGMEdge_ptr *d_edges = new RGMEdge_ptr[dsk*dsj*dsi*8];

  for(k=1,ik=0;k<dsk;k++,ik++){
    for(j=1,ij=0;j<dsj;j++,ij++){
      for(i=1,ii=0;i<dsi;i++,ii++){
        BCC_edges[ik*dsi*dsj*3 + ij*dsi*3 + ii*3 + 0] = 
          InsertEdge(
            BCC[(k-1)*(dsi+1)*(dsj+1) + j*(dsi+1) + i], 
            BCC[k*(dsi+1)*(dsj+1) + j*(dsi+1) + i]);
        BCC_edges[ik*dsi*dsj*3 + ij*dsi*3 + ii*3 + 1] = 
          InsertEdge(
            BCC[k*(dsi+1)*(dsj+1) + j*(dsi+1) + (i-1)], 
            BCC[k*(dsi+1)*(dsj+1) + j*(dsi+1) + i]);
        BCC_edges[ik*dsi*dsj*3 + ij*dsi*3 + ii*3 + 2] = 
          InsertEdge(
            BCC[k*(dsi+1)*(dsj+1) + (j-1)*(dsi+1) + i], 
            BCC[k*(dsi+1)*(dsj+1) + j*(dsi+1) + i]);

        // iBCC_edge[i] is parallel to BCC_edge[i]
        iBCC_edges[ik*dsi*dsj*3 + ij*dsi*3 + ii*3 + 0] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            iBCC[(ik+1)*dsi*dsj + ij*dsi + ii]);
        iBCC_edges[ik*dsi*dsj*3 + ij*dsi*3 + ii*3 + 1] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            iBCC[ik*dsi*dsj + ij*dsi + (ii+1)]);
        iBCC_edges[ik*dsi*dsj*3 + ij*dsi*3 + ii*3 + 2] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            iBCC[ik*dsi*dsj + (ij+1)*dsi + ii]);

        // reverse Z traversal of the BCC vertices above the current iBCC node
        d_edges[ik*dsi*dsj*8 + ij*dsi*8 + ii*8 + 0] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            BCC[k*(dsi+1)*(dsj+1) + j*(dsi+1) + i]);
        d_edges[ik*dsi*dsj*8 + ij*dsi*8 + ii*8 + 1] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            BCC[k*(dsi+1)*(dsj+1) + j*(dsi+1) + (i-1)]);
        d_edges[ik*dsi*dsj*8 + ij*dsi*8 + ii*8 + 2] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            BCC[k*(dsi+1)*(dsj+1) + (j-1)*(dsi+1) + i]);
        d_edges[ik*dsi*dsj*8 + ij*dsi*8 + ii*8 + 3] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            BCC[k*(dsi+1)*(dsj+1) + (j-1)*(dsi+1) + (i-1)]);


        // reverse Z traversal of the BCC vertices below the current iBCC node
        d_edges[ik*dsi*dsj*8 + ij*dsi*8 + ii*8 + 4] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            BCC[(k-1)*(dsi+1)*(dsj+1) + j*(dsi+1) + i]);
        d_edges[ik*dsi*dsj*8 + ij*dsi*8 + ii*8 + 5] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            BCC[(k-1)*(dsi+1)*(dsj+1) + j*(dsi+1) + (i-1)]);
        d_edges[ik*dsi*dsj*8 + ij*dsi*8 + ii*8 + 6] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            BCC[(k-1)*(dsi+1)*(dsj+1) + (j-1)*(dsi+1) + i]);
        d_edges[ik*dsi*dsj*8 + ij*dsi*8 + ii*8 + 7] = 
          InsertEdge(
            iBCC[ik*dsi*dsj + ij*dsi + ii], 
            BCC[(k-1)*(dsi+1)*(dsj+1) + (j-1)*(dsi+1) + (i-1)]);
      }
    }
  }

  for(k=1;k<dsk-2;k++){
    //      std::cout << "k=" << k << "(" << dsk-2 << ")" << std::endl;
    for(j=1;j<dsj-2;j++){
      //        std::cout << "j=" << k << "(" << dsj-2 << ")" << std::endl;
      for(i=1;i<dsi-2;i++){
        //          std::cout << "i=" << k << "(" << dsi-2 << ")" << std::endl;
        // 4 tetrahedra each sharing one edge of the 
        // iBCC cube
        RGMTetra_ptr new_tetra;

        // iBCC_edges[0]
        // edge 0 is opposite to 1, 
        // 2 -- to 3, and 4 -- to 5
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 1],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 0],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 0],   
          d_edges[(k+1)*dsi*dsj*8 + j*dsi*8 + i*8 + 5],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 1],   
          d_edges[(k+1)*dsi*dsj*8 + j*dsi*8 + i*8 + 4]);   
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 2],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 0],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 0],
          d_edges[(k+1)*dsi*dsj*8 + j*dsi*8 + i*8 + 6],
          d_edges[(k+1)*dsi*dsj*8 + j*dsi*8 + i*8 + 4],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 2]);
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + (j-1)*dsi*3 + i*3 + 1],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 0],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 2],
          d_edges[(k+1)*dsi*dsj*8 + j*dsi*8 + i*8 + 7],
          d_edges[(k+1)*dsi*dsj*8 + j*dsi*8 + i*8 + 6],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 3]);

        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + j*dsi*3 + (i-1)*3 + 2],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 0],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 1],
          d_edges[(k+1)*dsi*dsj*8 + j*dsi*8 + i*8 + 7],
          d_edges[(k+1)*dsi*dsj*8 + j*dsi*8 + i*8 + 5],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 3]);

        // iBCC_edges[1]
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 0],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 1],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 0],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + (i+1)*8 + 5],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + (i+1)*8 + 1],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 4]);
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 2],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 1],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 0],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + (i+1)*8 + 3],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + (i+1)*8 + 1],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 2]);
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + (j-1)*dsi*3 + i*3 + 0],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 1],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 2],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + (i+1)*8 + 7],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + (i+1)*8 + 3],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 6]);
        new_tetra = InsertTetraOriented(
          BCC_edges[(k-1)*dsi*dsj*3 + j*dsi*3 + i*3 + 2],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 1],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 4],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + (i+1)*8 + 7],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + (i+1)*8 + 5],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 6]);

        // iBCC_edges[2]
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 0],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 2],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 0],   
          d_edges[k*dsi*dsj*8 + (j+1)*dsi*8 + i*8 + 6],
          d_edges[k*dsi*dsj*8 + (j+1)*dsi*8 + i*8 + 2],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 4]);
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 1],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 2],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 0],   
          d_edges[k*dsi*dsj*8 + (j+1)*dsi*8 + i*8 + 3],
          d_edges[k*dsi*dsj*8 + (j+1)*dsi*8 + i*8 + 2],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 1]);
        new_tetra = InsertTetraOriented(
          BCC_edges[k*dsi*dsj*3 + j*dsi*3 + (i-1)*3 + 0],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 2],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 1],   
          d_edges[k*dsi*dsj*8 + (j+1)*dsi*8 + i*8 + 7],
          d_edges[k*dsi*dsj*8 + (j+1)*dsi*8 + i*8 + 3],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 5]);
        new_tetra = InsertTetraOriented(
          BCC_edges[(k-1)*dsi*dsj*3 + j*dsi*3 + i*3 + 1],
          iBCC_edges[k*dsi*dsj*3 + j*dsi*3 + i*3 + 2],
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 4],   
          d_edges[k*dsi*dsj*8 + (j+1)*dsi*8 + i*8 + 7],
          d_edges[k*dsi*dsj*8 + (j+1)*dsi*8 + i*8 + 6],   
          d_edges[k*dsi*dsj*8 + j*dsi*8 + i*8 + 5]);
      }
    }
  }

  delete [] BCC;
  delete [] iBCC;
  delete [] BCC_edges;
  delete [] iBCC_edges;
  delete [] d_edges;
}


/* Initial BCC spacing together with the number of resolutions ultimately
 * identify the final size of the mesh. It is also important to correctly
 * find the BCC spacing in order to represent all features of the mesh. */
template<class TInputImage, class TOutputMesh>
double
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::FindBCCSpacing(){
  double spacing;
  spacing = std::min(m_dimX, m_dimY);
  spacing = std::min(spacing, m_dimZ);
  return spacing/m_BCCSpacing;
}

/* Inserts a vertex */
template<class TInputImage, class TOutputMesh>
typename BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>::RGMVertex_ptr
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::InsertVertex(double cx, double cy, double cz){
  RGMVertex_ptr new_vertex = new RGMVertex;
  new_vertex->coords[0] = cx;
  new_vertex->coords[1] = cy;
  new_vertex->coords[2] = cz;
  m_Vertices.push_back(new_vertex);
  return new_vertex;
}

/* Inserts an edge */
template<class TInputImage, class TOutputMesh>
typename BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>::RGMEdge_ptr
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::InsertEdge(RGMVertex_ptr v0, RGMVertex_ptr v1){
  RGMEdge_ptr new_edge = new RGMEdge;
  new_edge->nodes[0] = v0;  
  new_edge->nodes[1] = v1;
  new_edge->children[0] = NULL;
  new_edge->children[0] = NULL;
  new_edge->midv = NULL;
  m_Edges.push_back(new_edge);
  return new_edge;
}

/* Insert tetrahedron, do the orientation check */
/* WARNING: only for the 0th level BCC!!! */
template<class TInputImage, class TOutputMesh>
typename BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>::RGMTetra_ptr
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::InsertTetraOriented(RGMEdge_ptr e0, RGMEdge_ptr e1, RGMEdge_ptr e2,
  RGMEdge_ptr e3, RGMEdge_ptr e4, RGMEdge_ptr e5){
  // Take e0 edge, edge e1 is opposite to e0.
  // Convention: imagine e0->nodes[0] as the front-most node of the
  // tetrahedron, then we want to order two edges sharing this vertex in
  // clockwise direction, identify their orientation (positive if
  // e_i->nodes[0] == e0->nodes[0]. For the edges of the face opposite to
  // e0->nodes[0] want to have them in clockwise direction too, same for
  // orientation.
  RGMTetra_ptr new_tetra = new RGMTetra;

  new_tetra->parent = NULL;
  new_tetra->subdiv = 0;
  new_tetra->edges_split = 0;
  new_tetra->edges_orient = 0;
  new_tetra->level = 0;
  new_tetra->greens = NULL;
  new_tetra->num_greens = 0;

  new_tetra->edges[0] = e0;
  e0->neihood.push_back(new_tetra);
  new_tetra->edges[1] = e1;
  e1->neihood.push_back(new_tetra);
  new_tetra->edges[2] = e2;
  e2->neihood.push_back(new_tetra);
  new_tetra->edges[3] = e3;
  e3->neihood.push_back(new_tetra);
  new_tetra->edges[4] = e4;
  e4->neihood.push_back(new_tetra);
  new_tetra->edges[5] = e5;
  e5->neihood.push_back(new_tetra);

  if(e2->nodes[0] != e0->nodes[0] &&
    e2->nodes[1] != e0->nodes[0]){
    // swap e2 and e3
    new_tetra->edges[2] = e3;
    new_tetra->edges[3] = e2;
  }
  if(e4->nodes[0] != e0->nodes[0] &&
    e4->nodes[1] != e0->nodes[0]){
    new_tetra->edges[4] = e5;
    new_tetra->edges[5] = e4;
  }

  // next we need to know the orientation of each of the edges relative to e0
  if(orient3d(&e0->nodes[0]->coords[0], &e0->nodes[1]->coords[0], 
      &e1->nodes[0]->coords[0], &e1->nodes[1]->coords[0])>0){ // DEST(e0)+e1 are in CW
    // e2 should share ORIG(e1)
    if(new_tetra->edges[4]->nodes[0] != e1->nodes[1] &&
      new_tetra->edges[4]->nodes[1] != e1->nodes[1]){
      // swap second with third pair of tetrahedron edges
      RGMEdge_ptr e_tmp;
      e_tmp = new_tetra->edges[2];
      new_tetra->edges[2] = new_tetra->edges[4];
      new_tetra->edges[4] = e_tmp;
      e_tmp = new_tetra->edges[3];
      new_tetra->edges[3] = new_tetra->edges[5];
      new_tetra->edges[5] = e_tmp;
    }
    if(new_tetra->edges[1]->nodes[1] != new_tetra->edges[3]->nodes[0])
      new_tetra->edges_orient |= (1<<3);
    if(new_tetra->edges[1]->nodes[1] != new_tetra->edges[5]->nodes[1])
      new_tetra->edges_orient |= (1<<5);
  } else {
    //std::cout << "orient3d()<0" << std::endl;
    new_tetra->edges_orient |= 2;
    // e2 should share DEST(e1)
    if(new_tetra->edges[4]->nodes[0] != e1->nodes[0] &&
      new_tetra->edges[4]->nodes[1] != e1->nodes[0]){
      // swap second with third pair of tetrahedron edges
      RGMEdge_ptr e_tmp;;
      e_tmp = new_tetra->edges[2];
      new_tetra->edges[2] = new_tetra->edges[4];
      new_tetra->edges[4] = e_tmp;
      e_tmp = new_tetra->edges[3];
      new_tetra->edges[3] = new_tetra->edges[5];
      new_tetra->edges[5] = e_tmp;
    }
    if(new_tetra->edges[1]->nodes[0] != new_tetra->edges[3]->nodes[0])
      new_tetra->edges_orient |= (1<<3);
    if(new_tetra->edges[1]->nodes[1] != new_tetra->edges[5]->nodes[1])
      new_tetra->edges_orient |= (1<<5);
  }
  // Now all edges of new_tetra are ordered:
  // e0, e2, and e4 share ORIG(e0) and go in CW,
  // and e1, e3, and e5 are in CW order. 

  // Find the orientation of each of the edges.
  // Edges sharing a node with e0 and have same 
  // ORIG that e0 have orientation 0.
  // Edges of the face opposite to ORIG(e0) and go in CW
  // have orientation 0. 
  // Otherwise orientation is 1.
  if(new_tetra->edges[0]->nodes[0] != new_tetra->edges[2]->nodes[0])
    new_tetra->edges_orient |= (1<<2);
  if(new_tetra->edges[0]->nodes[0] != new_tetra->edges[4]->nodes[0])
    new_tetra->edges_orient |= (1<<4);

  new_tetra->parent = NULL;
  new_tetra->subdiv = 0;
  new_tetra->edges_split = 0;

#ifdef VERBOSE
  std::cout << "Tetra #" << new_tetra->id << ": " <<
    new_tetra->edges[0] << "," << new_tetra->edges[1] << "," <<
    new_tetra->edges[2] << "," << new_tetra->edges[3] << "," <<
    new_tetra->edges[4] << "," << new_tetra->edges[5] << std::endl;
#endif

//  if(res>=m_current_level)
    m_Tetras.push_back(new_tetra);
//  else
//    m_pending_tetras.push_back(new_tetra);

  return new_tetra;
}

/* Go through the list of tetrahedra, and remove those which are completely
 * outside of the object's surface */
template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::RemoveBCCTetra(RGMTetra_ptr thisT){
  // This operation is safe, because tetrahedra do not have any descendants
  // yet. We go trhough all the edges and remove thisT from the
  // neighbourhoods. Edges with 0 neighbourhood lists can be deallocated.
  // Search is O(c) operation, because it's possible to bound the maximum
  // number of tets incident on an edge.
  for(int i=0;i<6;i++){
    typename std::list<RGMTetra_ptr>::iterator enI =
      std::find(thisT->edges[i]->neihood.begin(),
        thisT->edges[i]->neihood.end(), thisT);
    thisT->edges[i]->neihood.erase(enI);
    if(!thisT->edges[i]->neihood.size())
      delete thisT->edges[i];
  }
  delete thisT;
}

/* Return true if the tetrahedron has to be split based on the user specified
 * criteria */
template<class TInputImage, class TOutputMesh>
bool
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::TetraUnderrefined(RGMTetra_ptr thisT){
  if(m_SubdivInclusion)
    if(SubdivInclusionTest_Resolution(thisT, m_ResampleResolution))
      return true;
 
  for(typename std::list<SubdivisionTestFunctionPointer>::iterator
    scI=m_SubdivisionCriteria.begin();
    scI!=m_SubdivisionCriteria.end();
    scI++)
    if((*scI)(thisT->edges[0]->nodes[0]->coords, 
        thisT->edges[0]->nodes[1]->coords,
        thisT->edges[1]->nodes[0]->coords, 
        thisT->edges[1]->nodes[1]->coords, this))
      return true;
  
  return false;
}

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::FinalizeRedTetra(RGMTetra_ptr thisT){
  // 12 new face edges
  int i, face_id;
  RGMEdge_ptr  face_edges[4][3];
  RGMEdge_ptr  internal_edge;

  thisT->edges_split = 0x3F;

  // mark the tetrahedron and split each of its edges
  for(i=0;i<6;i++){

    if(!thisT->edges[i]->midv){
      SplitEdge(thisT->edges[i]);
    } else{
      assert(thisT->edges[i]->children[0]);
      assert(thisT->edges[i]->children[1]);
    }

    if(!thisT->edges[i]->children[0]){
      std::cout << "FAILURE, i=" << i << std::endl;
      assert(0);
    }
    assert(thisT->edges[i]->children[1]);
  } // for(all edges of thisT)

  if(!thisT->edges[1]->children[0]){
    std::cout << "FAILURE, i=" << i << std::endl;
    assert(0);
  }

  // choose the internal edge so that it is the shortest out of three
  // possible choices. This will create new tetras which are exactly those
  // as BCC tetras resulting from a mesh with cells 1/2 the size! [Molino]
  double internal_edge_length[3];
  internal_edge_length[0] = DistanceBwPoints(thisT->edges[0]->midv->coords, thisT->edges[1]->midv->coords);
  internal_edge_length[1] = DistanceBwPoints(thisT->edges[2]->midv->coords, thisT->edges[3]->midv->coords);
  internal_edge_length[2] = DistanceBwPoints(thisT->edges[4]->midv->coords, thisT->edges[5]->midv->coords);
  if(internal_edge_length[0]>internal_edge_length[1] &&
    internal_edge_length[0]>internal_edge_length[2]){
    RGMTetra rT;
    if(internal_edge_length[1]<internal_edge_length[2]){
      // reorder so that pair 2-3 is first
      rT.edges[0] = thisT->edges[2];
      rT.edges[1] = thisT->edges[3];
      rT.edges[2] = thisT->edges[4];
      rT.edges[3] = thisT->edges[5];
      rT.edges[4] = thisT->edges[0];
      rT.edges[5] = thisT->edges[1];
      rT.edges_orient = GetTetraEdgeOrient(&rT);
    } else {
      // reorder so that pair 4-5 is first
      rT.edges[0] = thisT->edges[4];
      rT.edges[1] = thisT->edges[5];
      rT.edges[2] = thisT->edges[0];
      rT.edges[3] = thisT->edges[1];
      rT.edges[4] = thisT->edges[2];
      rT.edges[5] = thisT->edges[3];
      rT.edges_orient = GetTetraEdgeOrient(&rT);
    }
    /*
    std::cout << "Need to reorder edges..." << std::endl;
    std::cout << "Lengths: " << internal_edge_length[0] << " ";
    std::cout << internal_edge_length[1] << " " << internal_edge_length[2] << std::endl;
    std::cout << "Tetra before: ";
     */
    thisT->edges[0] = rT.edges[0];
    thisT->edges[1] = rT.edges[1];
    thisT->edges[2] = rT.edges[2];
    thisT->edges[3] = rT.edges[3];
    thisT->edges[4] = rT.edges[4];
    thisT->edges[5] = rT.edges[5];
    thisT->edges_orient = rT.edges_orient;
  }
  internal_edge = InsertEdge(thisT->edges[0]->midv, thisT->edges[1]->midv);

  // 2) create 3 new edges for each of the faces, and put
  // them in a hashtable (order vertices of an edge new by 
  // pointer to facilitate lookup)
  // go through the faces
  for(face_id=0;face_id<4;face_id++){
    std::pair<RGMVertex_ptr,RGMVertex_ptr>  new_face_edge;
    unsigned char face_midv[3];
    // face_edges_order keeps all face edges listed in CW direction,
    // with the first edge opposing ORIG(e0) on the faces incident to 
    // ORIG(e0), and with the first edge parallel to e1 on the face opposite
    // to ORIG(e0)
    unsigned char face_edges_order[3] = {0,1,2};
    unsigned char tmp_midv;

    face_midv[0] = m_FaceToEdgesLT[face_id*3+0];
    face_midv[1] = m_FaceToEdgesLT[face_id*3+1];
    face_midv[2] = m_FaceToEdgesLT[face_id*3+2];

    // sort the face edges by the address of the midsection point
    if(thisT->edges[face_midv[0]]->midv <
      thisT->edges[face_midv[1]]->midv){   
      if(thisT->edges[face_midv[0]]->midv <
        thisT->edges[face_midv[2]]->midv){ 
        if(thisT->edges[face_midv[1]]->midv <
          thisT->edges[face_midv[2]]->midv){
          // do nothing
          // order of midpoints is {e0m, e2m, e5m}  {e4m, e0m, e3m}
          // edges are {e0m->e2m, e0m->e5m, e2m->e5m} {e4m->e0m, e4m->e3m, e0m->e3m}
          // order of face edges is {0,2,1}
          face_edges_order[1] = 2;
          face_edges_order[2] = 1;
        } else {  // b>=c
          // order of midpoints is {e0m, e5m, e2m} {e4m,e3m,e0m}
          // edges are {e0m->e5m, e0m->e2m, e5m->e2m} {e4m->e3m,e4m->e0m,e3m->e0m}
          // order of edges is {2,0,1}
          tmp_midv = face_midv[2];
          face_midv[2] = face_midv[1];
          face_midv[1] = tmp_midv;
          face_edges_order[0] = 2;
          face_edges_order[1] = 0;
          face_edges_order[2] = 1;
        }
      } else { // a>=c
        // order of midpoints is {e5m, e0m, e2m}
        // edges are {e5m->e0m, e5m->e2m, e0m->e2m}
        // order of edges is {2,1,0}
        tmp_midv = face_midv[0];
        face_midv[0] = face_midv[2];
        face_midv[2] = face_midv[1];
        face_midv[1] = tmp_midv;
        face_edges_order[0] = 2;
        face_edges_order[1] = 1;
        face_edges_order[2] = 0;
      }
    } else { // a>=b
      if(thisT->edges[face_midv[1]]->midv <
        thisT->edges[face_midv[2]]->midv){
        if(thisT->edges[face_midv[0]]->midv <
          thisT->edges[face_midv[2]]->midv){
          // order of midpoints is {e2m, e0m, e5m}
          // edges are {e2m->e0m, e2m->e5m, e0m->e5m}
          // order of edges is {0,1,2}
          tmp_midv = face_midv[0];
          face_midv[0] = face_midv[1];
          face_midv[1] = tmp_midv;
          face_edges_order[0] = 0;
          face_edges_order[1] = 1;
          face_edges_order[2] = 2;            
        } else { // a>=c
          // order of midpoints is {e2m, e5m, e0m}
          // edges are {e2m->e5m, e2m->e0m, e5m->e0m}
          // order is {1,0,2}
          tmp_midv = face_midv[0];
          face_midv[0] = face_midv[1];
          face_midv[1] = face_midv[2];
          face_midv[2] = tmp_midv;
          face_edges_order[0] = 1;
          face_edges_order[1] = 0;
          face_edges_order[2] = 2;            
        }
      } else { // b>=c
        // order of midpoints is {e5m, e2m, e0m} {e3m,e0m,e4m}
        // edges are {e5m->e2m, e5m->e0m, e2m->e0m} {e3m->e0m,e3m->e4m,e0m->e4m}
        // order is {2,1,0} {1,2,0}
        tmp_midv = face_midv[0];
        face_midv[0] = face_midv[2];
        face_midv[2] = tmp_midv;
        face_edges_order[0] = 1;    // this means, 0th edge which is 0->1 goes to the 2nd position
        face_edges_order[1] = 2;
        face_edges_order[2] = 0;            
      }
    }

    if(!thisT->edges[1]->children[0]){
      std::cout << "FAILURE! face_id=" << face_id << std::endl;
      assert(0);
    }

    face_edges[face_id][face_edges_order[0]] = 
      GetPutTmpEdge(thisT->edges[face_midv[0]]->midv, 
        thisT->edges[face_midv[1]]->midv);
    face_edges[face_id][face_edges_order[1]] = 
      GetPutTmpEdge(thisT->edges[face_midv[0]]->midv, 
        thisT->edges[face_midv[2]]->midv);
    face_edges[face_id][face_edges_order[2]] = 
      GetPutTmpEdge(thisT->edges[face_midv[1]]->midv, 
        thisT->edges[face_midv[2]]->midv);
  }

  // Finally, create internal tetrahedra...
  bool edge_orientations[6];
  for(i=0;i<6;i++)
    edge_orientations[i] = (thisT->edges_orient >> i) & 1;

  // By convention, children edges preserve parent orientation, and 
  // the first child is (ORIG(e) --> e->midv);
  // we also know the orientation of parent edges
  // AND we also know the relative positions of the new
  // edges, so that we don't need to reshuffle the edges'
  // pairs for new tetrahedra

  RGMTetra_ptr new_tetra;

  // Child #0
  // ORIG(e0) + 3 midpoints of edges incident to it
  new_tetra = InsertTetra(thisT->level+1,
    thisT->edges[0]->children[0],                    // first child of e0
    face_edges[1][0],                                // first edge in second face
    thisT->edges[2]->children[edge_orientations[2]], // child of e1
    face_edges[2][0],
    thisT->edges[4]->children[edge_orientations[4]],
    face_edges[0][0]);

  // orientation of edges that have parents is inherited from the parents
  new_tetra->edges_orient = (((unsigned char)edge_orientations[2])<<2) | 
    (((unsigned char)edge_orientations[4])<<4);
  // orientations of other edges has to be identified
  if(face_edges[1][0]->nodes[0] != thisT->edges[2]->midv)
    new_tetra->edges_orient |= (1<<1);
  if(face_edges[2][0]->nodes[0] != thisT->edges[4]->midv)
    new_tetra->edges_orient |= (1<<3);
  if(face_edges[0][0]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<5);

  // Child #1
  // DEST(e0) + 3 midpoints of edges incident to it
  new_tetra = InsertTetra(thisT->level+1,
    thisT->edges[0]->children[1], face_edges[3][0],
    face_edges[0][2], thisT->edges[3]->children[!edge_orientations[3]],
    face_edges[2][1], thisT->edges[5]->children[edge_orientations[5]]);

  new_tetra->edges_orient = (((unsigned char)edge_orientations[3])<<3) | 
    (((unsigned char)edge_orientations[5])<<5);
  if(face_edges[3][0]->nodes[0] != thisT->edges[5]->midv)
    new_tetra->edges_orient |= (1<<1);
  if(face_edges[0][2]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<2);
  if(face_edges[2][1]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<4);

  // Child #2
  // DEST(e2) + 3 midpoints of edges incident to it
  assert(thisT->edges[1]->children[0]);
  assert(thisT->edges[1]->children[1]);
  new_tetra = InsertTetra(thisT->level+1,
    face_edges[0][1], thisT->edges[1]->children[edge_orientations[1]],
    thisT->edges[2]->children[!edge_orientations[2]], face_edges[3][1],
    face_edges[1][2], thisT->edges[5]->children[!edge_orientations[5]]);

  if(new_tetra->edges[0]->nodes[0] == thisT->edges[2]->midv){
    new_tetra->edges_orient = 
      (((unsigned char)edge_orientations[1])<<1) | 
      (((unsigned char)edge_orientations[2])<<2) |
      (((unsigned char)edge_orientations[5])<<5);

    if(new_tetra->edges[3]->nodes[0] != thisT->edges[1]->midv)
      new_tetra->edges_orient |= (1<<3);
    if(new_tetra->edges[4]->nodes[0] != thisT->edges[2]->midv)
      new_tetra->edges_orient |= (1<<4);
  } else {
    RGMEdge_ptr tmp_edge;

    tmp_edge = new_tetra->edges[2];
    new_tetra->edges[2] = new_tetra->edges[3];
    new_tetra->edges[3] = tmp_edge;
    tmp_edge = new_tetra->edges[4];
    new_tetra->edges[4] = new_tetra->edges[5];
    new_tetra->edges[5] = tmp_edge;

    if(new_tetra->edges[1]->nodes[0] != thisT->edges[1]->midv)
      new_tetra->edges_orient |= (1<<1);
    if(new_tetra->edges[2]->nodes[0] != thisT->edges[5]->midv)
      new_tetra->edges_orient |= (1<<2);
    if(new_tetra->edges[3]->nodes[1] != thisT->edges[2]->midv)
      new_tetra->edges_orient |= (1<<3);
    if(new_tetra->edges[4]->nodes[0] != thisT->edges[5]->midv)
      new_tetra->edges_orient |= (1<<4);
    if(new_tetra->edges[5]->nodes[0] != thisT->edges[2]->midv)
      new_tetra->edges_orient |= (1<<5);
  }

  // Child #3
  // DEST(e4) + 3 midpoints of edges incident to it
  new_tetra = InsertTetra(thisT->level+1,
    face_edges[2][2], thisT->edges[1]->children[!edge_orientations[1]],
    face_edges[1][1], thisT->edges[3]->children[edge_orientations[3]],
    thisT->edges[4]->children[!edge_orientations[4]], face_edges[3][2]);

  if(new_tetra->edges[0]->nodes[0] == thisT->edges[4]->midv){
    new_tetra->edges_orient = (((unsigned char)edge_orientations[1])<<1) | 
      (((unsigned char)edge_orientations[3])<<3) |
      (((unsigned char)edge_orientations[4])<<4);

    if(new_tetra->edges[2]->nodes[0] != thisT->edges[4]->midv)
      new_tetra->edges_orient |= (1<<2);
    if(new_tetra->edges[5]->nodes[0] != thisT->edges[3]->midv)
      new_tetra->edges_orient |= (1<<5);
  } else {
    RGMEdge_ptr tmp_edge;

    tmp_edge = new_tetra->edges[2];
    new_tetra->edges[2] = new_tetra->edges[3];
    new_tetra->edges[3] = tmp_edge;
    tmp_edge = new_tetra->edges[4];
    new_tetra->edges[4] = new_tetra->edges[5];
    new_tetra->edges[5] = tmp_edge;

    if(new_tetra->edges[1]->nodes[1] != thisT->edges[1]->midv)
      new_tetra->edges_orient |= (1<<1);
    if(new_tetra->edges[2]->nodes[0] != thisT->edges[3]->midv)
      new_tetra->edges_orient |= (1<<2);
    if(new_tetra->edges[3]->nodes[0] != thisT->edges[1]->midv)
      new_tetra->edges_orient |= (1<<3);
    if(new_tetra->edges[4]->nodes[0] != thisT->edges[3]->midv)
      new_tetra->edges_orient |= (1<<4);
    if(new_tetra->edges[5]->nodes[0] != thisT->edges[4]->midv)
      new_tetra->edges_orient |= (1<<5);
  }

  // Child #4
  // internal_edge + midpoints of e2 and e4
  new_tetra = InsertTetra(thisT->level+1, 
    internal_edge, face_edges[1][0],
    face_edges[0][0], face_edges[1][1],
    face_edges[2][0], face_edges[1][2]);
  new_tetra->level = thisT->level+1;

  if(face_edges[1][0]->nodes[0] != thisT->edges[2]->midv)
    new_tetra->edges_orient |= (1<<1);
  if(face_edges[0][0]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<2);
  if(face_edges[1][1]->nodes[0] != thisT->edges[4]->midv)
    new_tetra->edges_orient |= (1<<3);
  if(face_edges[2][0]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<4);
  if(face_edges[1][2]->nodes[0] != thisT->edges[1]->midv)
    new_tetra->edges_orient |= (1<<5);

  // Child #5
  // internal_edge + midpoints of e4 and e3
  new_tetra = InsertTetra(thisT->level+1,
    internal_edge, face_edges[2][2],
    face_edges[2][0], face_edges[3][2],
    face_edges[2][1], face_edges[1][1]);
  new_tetra->level = thisT->level+1;

  if(face_edges[2][2]->nodes[0] != thisT->edges[4]->midv)
    new_tetra->edges_orient |= (1<<1);
  if(face_edges[2][0]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<2);
  if(face_edges[3][2]->nodes[0] != thisT->edges[3]->midv)
    new_tetra->edges_orient |= (1<<3);
  if(face_edges[2][1]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<4);
  if(face_edges[1][1]->nodes[0] != thisT->edges[1]->midv)
    new_tetra->edges_orient |= (1<<5);

  // Child #6
  // internal_edge + midpoints of e3 and e5
  new_tetra = InsertTetra(thisT->level+1,
    internal_edge, face_edges[3][0],
    face_edges[2][1], face_edges[3][1],
    face_edges[0][2], face_edges[3][2]);
  new_tetra->level = thisT->level+1;

  if(face_edges[3][0]->nodes[0] != thisT->edges[3]->midv)
    new_tetra->edges_orient |= (1<<1);
  if(face_edges[2][1]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<2);
  if(face_edges[3][1]->nodes[0] != thisT->edges[5]->midv)
    new_tetra->edges_orient |= (1<<3);
  if(face_edges[0][2]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<4);
  if(face_edges[3][2]->nodes[0] != thisT->edges[1]->midv)
    new_tetra->edges_orient |= (1<<5);

  // Child #7
  // internal_edge + midpoints of e2 and e5
  new_tetra = InsertTetra(thisT->level+1,
    internal_edge, face_edges[0][1],
    face_edges[0][2], face_edges[1][2],
    face_edges[0][0], face_edges[3][1]);
  new_tetra->level = thisT->level+1;

  if(face_edges[0][1]->nodes[0] != thisT->edges[5]->midv)
    new_tetra->edges_orient |= (1<<1);
  if(face_edges[0][2]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<2);
  if(face_edges[1][2]->nodes[0] != thisT->edges[2]->midv)
    new_tetra->edges_orient |= (1<<3);
  if(face_edges[0][0]->nodes[0] != thisT->edges[0]->midv)
    new_tetra->edges_orient |= (1<<4);
  if(face_edges[3][1]->nodes[0] != thisT->edges[1]->midv)
    new_tetra->edges_orient |= (1<<5);

  RemoveTetra(thisT);
}

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::FinalizeGreenTetra(RGMTetra_ptr thisT){
  int split_edges_ids[3];
  int edges_split_total;
  int first_split_edge;
  int split_edge_orient, opposite_edge;
  int i, j;

  assert(thisT->greens==NULL);
  assert(thisT->parent==NULL);

  edges_split_total = 0;
  first_split_edge = -1;
  // identify the configuration
  for(i=0;i<6;i++){
    if(thisT->edges[i]->midv){
      split_edges_ids[edges_split_total++] = i;
      if(first_split_edge==-1)
        first_split_edge = i;
      thisT->edges_split |= 1<<i;
    }
  }
  //    assert(edges_split_total==TetraEdgeConf(thisT));

  opposite_edge = first_split_edge ^ 1;
  split_edge_orient = (thisT->edges_orient>>first_split_edge) & 1;

  switch(edges_split_total){
  case 0:{
    assert(0);
    break;
  }
  case 1:{
    //std::cout << "Green: 1" << std::endl;
    // need to create two face edges, which connect the midpoint with the
    // vertices of the opposite edge (unless they have already been
    // created). Flip the last bit to get the opposite edge.
    RGMEdge_ptr face_edges[2];
    RGMTetra_ptr new_tetra;

    thisT->greens = new RGMTetra[2];

      face_edges[0] = GetPutTmpEdge(thisT->edges[opposite_edge]->nodes[0],
        thisT->edges[first_split_edge]->midv);
      face_edges[1] = GetPutTmpEdge(thisT->edges[opposite_edge]->nodes[1],
        thisT->edges[first_split_edge]->midv);

      // disregard order of edges, because GREEN children will not be split
      new_tetra = InsertTetraAllocated(thisT->level+1,
        thisT->edges[first_split_edge]->children[0], thisT->edges[opposite_edge],
        thisT->edges[m_IncidentEdgesLT[12*split_edge_orient+first_split_edge*2]], face_edges[0],
        thisT->edges[m_IncidentEdgesLT[12*split_edge_orient+first_split_edge*2+1]], face_edges[1],
        &thisT->greens[0]);
      new_tetra->parent = thisT;

    new_tetra = InsertTetraAllocated(thisT->level+1,
      thisT->edges[first_split_edge]->children[1], thisT->edges[opposite_edge],
      face_edges[0], thisT->edges[m_IncidentEdgesLT[12*(!split_edge_orient)+first_split_edge*2]],
      face_edges[1], thisT->edges[m_IncidentEdgesLT[12*(!split_edge_orient)+first_split_edge*2+1]],
      &thisT->greens[1]);
    new_tetra->parent = thisT;

    thisT->num_greens = 2;
    break;
  }
  case 2:{
    //std::cout << "Green:2" << std::endl;
    RGMEdge_ptr face_edges[4];
    RGMEdge_ptr internal_edge;
    RGMTetra_ptr new_tetra;
    int opposite_edge = first_split_edge ^ 1;

    thisT->greens = new RGMTetra[4];

    if(split_edges_ids[0] != (split_edges_ids[1]^1))
      assert(0);

    assert(thisT->edges[first_split_edge]->children[0]);
    assert(thisT->edges[opposite_edge]->children[0]);

    face_edges[0] = GetPutTmpEdge(thisT->edges[opposite_edge]->nodes[0],
      thisT->edges[first_split_edge]->midv);
    face_edges[1] = GetPutTmpEdge(thisT->edges[opposite_edge]->nodes[1],
      thisT->edges[first_split_edge]->midv);
    face_edges[2] = GetPutTmpEdge(thisT->edges[first_split_edge]->nodes[0],
      thisT->edges[opposite_edge]->midv);
    face_edges[3] = GetPutTmpEdge(thisT->edges[first_split_edge]->nodes[1],
      thisT->edges[opposite_edge]->midv);

    internal_edge = InsertEdge(thisT->edges[first_split_edge]->midv,
      thisT->edges[opposite_edge]->midv);

    if( (((thisT->edges_orient)>>first_split_edge)&1) ^ // this is the orientation of opposite_edge 
      (((thisT->edges_orient)>>opposite_edge)&1) ){   //   relative to first_split_edge

      new_tetra = InsertTetraAllocated(thisT->level+1,
        internal_edge, thisT->edges[m_IncidentEdgesLT[12*split_edge_orient+first_split_edge*2]],
        face_edges[2], face_edges[1],
        thisT->edges[opposite_edge]->children[1], thisT->edges[first_split_edge]->children[0],
        &thisT->greens[0]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        internal_edge, thisT->edges[m_IncidentEdgesLT[12*split_edge_orient+first_split_edge*2+1]],
        thisT->edges[opposite_edge]->children[0], thisT->edges[first_split_edge]->children[0],
        face_edges[2], face_edges[0],
        &thisT->greens[1]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        internal_edge, thisT->edges[m_IncidentEdgesLT[12*(!split_edge_orient)+first_split_edge*2]],
        face_edges[3], face_edges[0],
        thisT->edges[opposite_edge]->children[0], thisT->edges[first_split_edge]->children[1],
        &thisT->greens[2]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        internal_edge, thisT->edges[m_IncidentEdgesLT[12*(!split_edge_orient)+first_split_edge*2+1]],
        thisT->edges[opposite_edge]->children[1], thisT->edges[first_split_edge]->children[1],
        face_edges[3], face_edges[1],
        &thisT->greens[3]);
      new_tetra->parent = thisT;
    } else {
      new_tetra = InsertTetraAllocated(thisT->level+1,
        internal_edge, thisT->edges[m_IncidentEdgesLT[12*split_edge_orient+first_split_edge*2]],
        face_edges[2], face_edges[0],
        thisT->edges[opposite_edge]->children[0], thisT->edges[first_split_edge]->children[0],
        &thisT->greens[0]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        internal_edge, thisT->edges[m_IncidentEdgesLT[12*split_edge_orient+first_split_edge*2+1]],
        thisT->edges[opposite_edge]->children[1], thisT->edges[first_split_edge]->children[0],
        face_edges[2], face_edges[1],
        &thisT->greens[1]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        internal_edge, thisT->edges[m_IncidentEdgesLT[12*(!split_edge_orient)+first_split_edge*2]],
        face_edges[3], face_edges[1],
        thisT->edges[opposite_edge]->children[1], thisT->edges[first_split_edge]->children[1],
        &thisT->greens[2]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        internal_edge, thisT->edges[m_IncidentEdgesLT[12*(!split_edge_orient)+first_split_edge*2+1]],
        thisT->edges[opposite_edge]->children[0], thisT->edges[first_split_edge]->children[1],
        face_edges[3], face_edges[0],
        &thisT->greens[3]);
      new_tetra->parent = thisT;
    }
    thisT->num_greens = 4;
    break;
  }
  case 3:{
    //std::cout << "Green: 3" << std::endl;
    RGMEdge_ptr long_face_edges[3];
    RGMEdge_ptr short_face_edges[3];
    RGMTetra_ptr new_tetra;
    RGMVertex_ptr opposite_vertex;

    assert(split_edges_ids[0] == 
      m_ThirdFaceEdgeLT[split_edges_ids[1]*6+split_edges_ids[2]]);

    thisT->greens = new RGMTetra[4];

    // which 3 edges are split?
    split_edges_ids[0] = first_split_edge;

    RGMTetra tmp_tetra;
    bool tmp_orientations[6];

    tmp_tetra.edges[0] = thisT->edges[split_edges_ids[0]];
    tmp_tetra.edges[1] = thisT->edges[split_edges_ids[0]^1];

    tmp_tetra.edges[2] = thisT->edges[m_IncidentEdgesLT[12*split_edge_orient+split_edges_ids[0]*2]];
    tmp_tetra.edges[3] = thisT->edges[m_IncidentEdgesLT[12*(!split_edge_orient)+split_edges_ids[0]*2]];
    tmp_tetra.edges[4] = thisT->edges[m_IncidentEdgesLT[12*split_edge_orient+split_edges_ids[0]*2+1]];
    tmp_tetra.edges[5] = thisT->edges[m_IncidentEdgesLT[12*(!split_edge_orient)+split_edges_ids[0]*2+1]];

    tmp_tetra.edges_split = 0;
    tmp_tetra.parent = NULL;

    tmp_orientations[1] = ((thisT->edges_orient>>split_edges_ids[0])&1) ^
      ((thisT->edges_orient>>split_edges_ids[1])&1);
    tmp_orientations[2] = !(tmp_tetra.edges[2]->nodes[0]==tmp_tetra.edges[0]->nodes[0]);
    tmp_orientations[3] = !(tmp_tetra.edges[3]->nodes[1]==tmp_tetra.edges[0]->nodes[1]);
    tmp_orientations[4] = !(tmp_tetra.edges[4]->nodes[0]==tmp_tetra.edges[0]->nodes[0]);
    tmp_orientations[5] = !(tmp_tetra.edges[5]->nodes[0]==tmp_tetra.edges[0]->nodes[1]);
    
    tmp_tetra.edges_orient = 0;
    for(j=0;j<6;j++)
      tmp_tetra.edges_orient |= (tmp_orientations[j]<<j);

    if(tmp_tetra.edges[2]->midv){

      assert(tmp_tetra.edges[2]->midv);
      assert(tmp_tetra.edges[5]->midv);

      opposite_vertex = tmp_tetra.edges[4]->nodes[!tmp_orientations[4]];

      long_face_edges[0] = GetPutTmpEdge(opposite_vertex, tmp_tetra.edges[0]->midv);
      long_face_edges[1] = GetPutTmpEdge(opposite_vertex, tmp_tetra.edges[2]->midv);
      long_face_edges[2] = GetPutTmpEdge(opposite_vertex, tmp_tetra.edges[5]->midv);

      // Should it be GetTmpEdge?..
      short_face_edges[0] = GetPutTmpEdge(tmp_tetra.edges[0]->midv, tmp_tetra.edges[2]->midv);
      short_face_edges[1] = GetPutTmpEdge(tmp_tetra.edges[0]->midv, tmp_tetra.edges[5]->midv);
      short_face_edges[2] = GetPutTmpEdge(tmp_tetra.edges[2]->midv, tmp_tetra.edges[5]->midv);

      new_tetra = InsertTetraAllocated(thisT->level+1,
        tmp_tetra.edges[1], short_face_edges[2],
        tmp_tetra.edges[2]->children[!tmp_orientations[2]], long_face_edges[2],
        tmp_tetra.edges[5]->children[!tmp_orientations[5]], long_face_edges[1],
        &thisT->greens[0]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        tmp_tetra.edges[3], short_face_edges[1],
        tmp_tetra.edges[5]->children[tmp_orientations[5]], long_face_edges[0],
        tmp_tetra.edges[0]->children[1], long_face_edges[2],
        &thisT->greens[1]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        tmp_tetra.edges[4], short_face_edges[0],
        tmp_tetra.edges[0]->children[0], long_face_edges[1],
        tmp_tetra.edges[2]->children[tmp_orientations[2]], long_face_edges[0],
        &thisT->greens[2]);
      new_tetra->parent = thisT;
    } else {

      if(!tmp_tetra.edges[4]->midv){
        assert(0);
      }
      assert(tmp_tetra.edges[3]->midv);

      opposite_vertex = tmp_tetra.edges[2]->nodes[!tmp_orientations[2]];

      long_face_edges[0] = GetPutTmpEdge(opposite_vertex, tmp_tetra.edges[0]->midv);
      long_face_edges[1] = GetPutTmpEdge(opposite_vertex, tmp_tetra.edges[4]->midv);
      long_face_edges[2] = GetPutTmpEdge(opposite_vertex, tmp_tetra.edges[3]->midv);

      // Should it be GetTmpEdge?..
      short_face_edges[0] = GetPutTmpEdge(tmp_tetra.edges[0]->midv, tmp_tetra.edges[4]->midv);
      short_face_edges[1] = GetPutTmpEdge(tmp_tetra.edges[0]->midv, tmp_tetra.edges[3]->midv);
      short_face_edges[2] = GetPutTmpEdge(tmp_tetra.edges[4]->midv, tmp_tetra.edges[3]->midv);

      new_tetra = InsertTetraAllocated(thisT->level+1,
        tmp_tetra.edges[1], short_face_edges[2],
        tmp_tetra.edges[4]->children[!tmp_orientations[4]], long_face_edges[2],
        tmp_tetra.edges[3]->children[tmp_orientations[3]], long_face_edges[1],
        &thisT->greens[0]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1, 
        tmp_tetra.edges[5], short_face_edges[1],
        tmp_tetra.edges[3]->children[!tmp_orientations[3]], long_face_edges[0],
        tmp_tetra.edges[0]->children[1], long_face_edges[2],
        &thisT->greens[1]);
      new_tetra->parent = thisT;

      new_tetra = InsertTetraAllocated(thisT->level+1,
        tmp_tetra.edges[2], short_face_edges[0],
        tmp_tetra.edges[0]->children[0], long_face_edges[1],
        tmp_tetra.edges[4]->children[tmp_orientations[4]], long_face_edges[0],
        &thisT->greens[2]);
      new_tetra->parent = thisT;
    }
    new_tetra = InsertTetraAllocated(thisT->level+1,
      long_face_edges[0], short_face_edges[2],
      long_face_edges[1], short_face_edges[1],
      long_face_edges[2], short_face_edges[0],
      &thisT->greens[3]);
    new_tetra->parent = thisT;

    thisT->num_greens = 4;
    break;
  }
  default:{
    std::cout << edges_split_total << " edges split.... Tetra type is " << 
      ((int)thisT->edges_split) << std::endl;
    assert(0);
  }
  }

}

template<class TInputImage, class TOutputMesh>
typename BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>::RGMEdge_ptr
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::GetPutTmpEdge(RGMVertex_ptr v0, RGMVertex_ptr v1){
  RGMEdge_ptr face_edge = NULL;
  std::pair<RGMVertex_ptr, RGMVertex_ptr> v_pair;
  typename std::map<std::pair<RGMVertex_ptr, RGMVertex_ptr>, RGMEdge_ptr,ltVertexPair>::iterator edge_mapI;

  if(v0>v1){
    v_pair.first = v1;
    v_pair.second = v0;
  } else {
    v_pair.first = v0;
    v_pair.second = v1;
  }

  if((edge_mapI = m_TmpEdgeMap.find(v_pair)) == m_TmpEdgeMap.end()){
    face_edge = InsertEdge(v_pair.first, v_pair.second);
    m_TmpEdgeMap[v_pair] = face_edge;
  } else {
    face_edge = (*edge_mapI).second;
    if((face_edge->nodes[0] != v_pair.first) ||
      (face_edge->nodes[1] != v_pair.second)){
      std::cout << "Edge ptr: " << face_edge << std::endl;
      std::cout << "Lookup nodes: " << v_pair.first << ", " 
        << v_pair.second << std::endl;
      std::cout << "Edge nodes: " << face_edge->nodes[0] << ", " 
        << face_edge->nodes[1] << std::endl;
      std::cout << "Neighbors: ";
      for(typename std::list<RGMTetra_ptr>::iterator nlI = face_edge->neihood.begin();
        nlI!=face_edge->neihood.end();nlI++)
        std::cout << *nlI << " ";
      std::cout << std::endl;
      assert(0);
    }
    m_TmpEdgeMap.erase(edge_mapI);
  }
  return face_edge;
}

template<class TInputImage, class TOutputMesh>
typename BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>::RGMTetra_ptr
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::InsertTetra(unsigned char resolution, 
  RGMEdge_ptr e0, RGMEdge_ptr e1, RGMEdge_ptr e2, 
  RGMEdge_ptr e3, RGMEdge_ptr e4, RGMEdge_ptr e5){
  // assume that edge orientation is correct
  RGMTetra_ptr new_tetra = new RGMTetra;

  new_tetra->parent = NULL;
  new_tetra->subdiv = 0;
  new_tetra->edges_split = 0;
  new_tetra->edges_orient = 0;
  new_tetra->level = resolution;
  new_tetra->greens = NULL;
  new_tetra->num_greens = 0;

  assert(e0);
  assert(e1);
  assert(e2);
  assert(e3);
  assert(e4);
  assert(e5);

  new_tetra->edges[0] = e0;
  e0->neihood.push_back(new_tetra);
  new_tetra->edges[1] = e1;
  e1->neihood.push_back(new_tetra);
  new_tetra->edges[2] = e2;
  e2->neihood.push_back(new_tetra);
  new_tetra->edges[3] = e3;
  e3->neihood.push_back(new_tetra);
  new_tetra->edges[4] = e4;
  e4->neihood.push_back(new_tetra);
  new_tetra->edges[5] = e5;
  e5->neihood.push_back(new_tetra);

  new_tetra->parent = NULL;
  new_tetra->subdiv = 0;
  new_tetra->edges_split = 0;

  if(resolution>=m_CurrentResol)
    m_Tetras.push_back(new_tetra);
  else
    m_PendingTetras.push_back(new_tetra);

  return new_tetra;
}

template<class TInputImage, class TOutputMesh>
unsigned char
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::GetTetraEdgeConf(RGMTetra_ptr thisT){
  unsigned char edges_split_total = 0;

  thisT->edges_split = 0;
  for(int i=0;i<6;i++)
    if(thisT->edges[i]->midv){
      thisT->edges_split |= (1<<i);
      edges_split_total++;
    }
  return edges_split_total;
}

template<class TInputImage, class TOutputMesh>
unsigned char
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::GetTetraEdgeOrient(RGMTetra_ptr thisT){
  RGMEdge_ptr e0, e1, e2, e3, e4, e5;
  
  e0 = thisT->edges[0];
  e1 = thisT->edges[1];
  e2 = thisT->edges[2];
  e3 = thisT->edges[3];
  e4 = thisT->edges[4];
  e5 = thisT->edges[5];
  thisT->edges_orient = 0;   

  if(e2->nodes[0] != e0->nodes[0] &&
    e2->nodes[1] != e0->nodes[0]){
    // swap e2 and e3
    thisT->edges[2] = e3;
    thisT->edges[3] = e2;
  }
  if(e4->nodes[0] != e0->nodes[0] &&
    e4->nodes[1] != e0->nodes[0]){
    thisT->edges[4] = e5;
    thisT->edges[5] = e4;
  }

  if(orient3d(&(e0->nodes[0]->coords[0]), 
      &(e0->nodes[1]->coords[0]), 
      &(e1->nodes[0]->coords[0]), 
      &(e1->nodes[1]->coords[0])) > 0){ 

    if(thisT->edges[4]->nodes[0] != e1->nodes[1] &&
      thisT->edges[4]->nodes[1] != e1->nodes[1]){
      // swap second with third pair of tetrahedron edges
      RGMEdge_ptr e_tmp;
      e_tmp = thisT->edges[2];
      thisT->edges[2] = thisT->edges[4];
      thisT->edges[4] = e_tmp;
      e_tmp = thisT->edges[3];
      thisT->edges[3] = thisT->edges[5];
      thisT->edges[5] = e_tmp;
    }

    if(thisT->edges[1]->nodes[1] != thisT->edges[3]->nodes[0])
      thisT->edges_orient |= (1<<3);
    if(thisT->edges[1]->nodes[0] != thisT->edges[5]->nodes[1])
      thisT->edges_orient |= (1<<5);
  } else {
    thisT->edges_orient |= 2;
    if(thisT->edges[4]->nodes[0] != e1->nodes[0] &&
      thisT->edges[4]->nodes[1] != e1->nodes[0]){
      // swap second with third pair of tetrahedron edges
      RGMEdge_ptr e_tmp;;
      e_tmp = thisT->edges[2];
      thisT->edges[2] = thisT->edges[4];
      thisT->edges[4] = e_tmp;
      e_tmp = thisT->edges[3];
      thisT->edges[3] = thisT->edges[5];
      thisT->edges[5] = e_tmp;
    }

    if(thisT->edges[1]->nodes[0] != thisT->edges[3]->nodes[0])
      thisT->edges_orient |= (1<<3);
    if(thisT->edges[1]->nodes[1] != thisT->edges[5]->nodes[1])
      thisT->edges_orient |= (1<<5);
  }
  // Now all edges of new_tetra are ordered:
  // e0, e2, and e4 share ORIG(e0) and go in CW,
  // and e1, e3, and e5 are in CW order. 

  // Find the orientation of each of the edges.
  // Edges sharing a node with e0 and have same 
  // ORIG that e0 have orientation 0.
  // Edges of the face opposite to ORIG(e0) and go in CW
  // have orientation 0. Otherwise orientation is 1.
  if(thisT->edges[0]->nodes[0] != thisT->edges[2]->nodes[0])
    thisT->edges_orient |= (1<<2);
  if(thisT->edges[0]->nodes[0] != thisT->edges[4]->nodes[0])
    thisT->edges_orient |= (1<<4);
  return thisT->edges_orient;
}

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::RemoveTetraAllocated(RGMTetra_ptr thisT){
  int i;

  for(i=0;i<6;i++){
    typename std::list<RGMTetra_ptr>::iterator enI = std::find(thisT->edges[i]->neihood.begin(),
      thisT->edges[i]->neihood.end(), thisT);
    assert(enI!=thisT->edges[i]->neihood.end());
    thisT->edges[i]->neihood.erase(enI);
  }
}

template<class TInputImage, class TOutputMesh>
typename BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>::RGMTetra_ptr
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::InsertTetraAllocated(int res, 
  RGMEdge_ptr e0, RGMEdge_ptr e1, RGMEdge_ptr e2,
  RGMEdge_ptr e3, RGMEdge_ptr e4, RGMEdge_ptr e5,
  RGMTetra_ptr thisT){
  // assume that edge orientation is correct
  RGMTetra_ptr new_tetra = thisT;

  new_tetra->parent = NULL;
  new_tetra->subdiv = 0;
  new_tetra->edges_split = 0;
  new_tetra->edges_orient = 0;
  new_tetra->level = res;
  new_tetra->greens = NULL;
  new_tetra->num_greens = 0;

  assert(e0);
  assert(e1);
  assert(e2);
  assert(e3);
  assert(e4);
  assert(e5);

  new_tetra->edges[0] = e0;
  e0->neihood.push_back(new_tetra);
  new_tetra->edges[1] = e1;
  e1->neihood.push_back(new_tetra);
  new_tetra->edges[2] = e2;
  e2->neihood.push_back(new_tetra);
  new_tetra->edges[3] = e3;
  e3->neihood.push_back(new_tetra);
  new_tetra->edges[4] = e4;
  e4->neihood.push_back(new_tetra);
  new_tetra->edges[5] = e5;
  e5->neihood.push_back(new_tetra);

  new_tetra->parent = NULL;
  new_tetra->subdiv = 0;
  new_tetra->edges_split = 0;
  new_tetra->level = res;

  if(res>=(int)m_CurrentResol)
    m_Tetras.push_back(new_tetra);
  else
    m_PendingTetras.push_back(new_tetra);

  return new_tetra;
}

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::RemoveTetra(RGMTetra_ptr thisT){
  int i;
    
  for(i=0;i<6;i++){
    typename std::list<RGMTetra_ptr>::iterator enI = std::find(thisT->edges[i]->neihood.begin(),
      thisT->edges[i]->neihood.end(), thisT);
    assert(enI!=thisT->edges[i]->neihood.end());
    thisT->edges[i]->neihood.erase(enI);
  }
  delete thisT;
}

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::SplitEdge(RGMEdge_ptr thisE){
  RGMVertex_ptr newV = new RGMVertex;

  newV->coords[0] = (thisE->nodes[0]->coords[0] + 
    thisE->nodes[1]->coords[0])/2.0;
  newV->coords[1] = (thisE->nodes[0]->coords[1] +
    thisE->nodes[1]->coords[1])/2.0;
  newV->coords[2] = (thisE->nodes[0]->coords[2] +
    thisE->nodes[1]->coords[2])/2.0;
  thisE->midv = newV;
  m_Vertices.push_back(newV);

  thisE->children[0] = 
    InsertEdge(thisE->nodes[0], thisE->midv);
  thisE->children[1] = 
    InsertEdge(thisE->midv, thisE->nodes[1]);
}

template<class TInputImage, class TOutputMesh>
void
BinaryMaskTo3DAdaptiveMeshFilter<TInputImage,TOutputMesh>
::PrepareInputImage()
 {
  typename CastFilterType::Pointer cast_filter = CastFilterType::New();
  typename ResamplerType::Pointer resampler = ResamplerType::New();
  IdentityTransformType::Pointer transform = IdentityTransformType::New();
  typename InputImageType::SpacingType input_spacing = this->m_InputImage->GetSpacing();
  typename InputImageType::SpacingType output_spacing;
  typename InputImageType::SizeType input_size = this->m_InputImage->GetLargestPossibleRegion().GetSize();
  typename InputImageType::SizeType output_size;
  typedef NearestNeighborInterpolateImageFunction<InternalImageType, double> ResampleInterpolatorType;
  typename ResampleInterpolatorType::Pointer resample_interpolator = ResampleInterpolatorType::New();

  typename InputImageType::DirectionType identityDirection;
  identityDirection.SetIdentity();
  if(this->m_InputImage->GetDirection() != identityDirection)
 	  itkExceptionMacro("The direction of the input image is : " << this->m_InputImage->GetDirection() << " Should be :" << identityDirection);

  // 1) Cast the input to the internal image type
  cast_filter->SetInput(this->m_InputImage);
  cast_filter->Update();

  // 2) Resample the input image to unit sized voxels
  output_spacing[0] = 1.0;
  output_spacing[1] = 1.0;
  output_spacing[2] = 1.0;

  output_size[0] = static_cast<typename InputImageType::SizeType::SizeValueType>
    (ceil((double)input_size[0]*input_spacing[0]));
  output_size[1] = static_cast<typename InputImageType::SizeType::SizeValueType>
    (ceil((double)input_size[1]*input_spacing[1]));
  output_size[2] = static_cast<typename InputImageType::SizeType::SizeValueType>
    (ceil((double)input_size[2]*input_spacing[2]));

  transform->SetIdentity();
  resampler->SetTransform(transform);
  resampler->SetInterpolator(resample_interpolator);
  resampler->SetOutputSpacing(output_spacing);
  resampler->SetOutputOrigin(this->m_InputImage->GetOrigin());
  resampler->SetSize(output_size);
  resampler->SetInput(cast_filter->GetOutput());
  resampler->Update();

  this->m_ReadyInputImage = resampler->GetOutput();
  //std::cout << "Input image was resampled " << std::endl;
  //std::cout<<"Size: " << this->m_ReadyInputImage->GetLargestPossibleRegion().GetSize() << std::endl;
  //std::cout<<"Spacing: " << this->m_ReadyInputImage->GetSpacing() << std::endl;
  //std::cout<<"Origin: " << this->m_ReadyInputImage->GetOrigin() << std::endl;
  //std::cout<<"Direction: " << this->m_ReadyInputImage->GetDirection() << std::endl;
}


} /** end namespace itk. */

#endif
