/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "T2PatchConstructor.hpp"
#include "../Helper/HalfedgeOperation.hpp"
#include "../Helper/ReadCSV2Matrix.hpp"
#include "../Patch/Patch.hpp"

Mat256x20d T2PatchConstructor::getMask()
{
    std::string t_MaskCSVFilePathT2 = "./Table/T2.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathT2, 256, 20); //  NumOfTotalCpts = 256, NumOfVerts = 20
}


/*
 *    T2 Structure
 *
 *    0 - 1 -------- 2 - 3
 *    |   |          |   |
 *    4 - 5 -------- 6 - 7
 *    |   |          |   |
 *    8 - 9          |   |
 *    |   |          |   |
 *   10 - 11 - 12 - 13 - 14
 *    |   |    |     |   |
 *   15 - 16 - 17 - 18 - 19
 *
 */
bool T2PatchConstructor::isSamePatchType(const FaceHandle& a_FaceHandle)
{
    // Check if there are 6 vertices of face
    if(!Helper::is_hexagon(m_Mesh, a_FaceHandle))
    {
        return false;
    }

    // Get vert 13 6 5 9 11 12 and store them to vector
    auto t_VertsAroundFace = Helper::get_verts_of_face(m_Mesh, a_FaceHandle);

    // Get valence of surounding vertices
    std::vector<int> t_ValencesOfVerts;
    for(auto t_Vert : t_VertsAroundFace)
    {
        t_ValencesOfVerts.push_back(Helper::get_vert_valence(m_Mesh, t_Vert));
    }

    // Rotate vector untill the first element is 3
    for(int i=0; i<t_ValencesOfVerts.size(); i++)
    {
        const int t_ThreeValence = 3;
        if(t_ValencesOfVerts[0]==t_ThreeValence)
        {
            break;
        }
        std::rotate(t_ValencesOfVerts.begin(),
                    t_ValencesOfVerts.begin() + 1, t_ValencesOfVerts.end());
    }

    // After rotation the order of valence should be 3 4 4 4 3 4 or 3 4 3 4 4 4
    std::vector<int> t_CorrectOrder1 = {3,4,4,4,3,4};
    std::vector<int> t_CorrectOrder2 = {3,4,3,4,4,4};
    if(t_ValencesOfVerts!=t_CorrectOrder1 && t_ValencesOfVerts!=t_CorrectOrder2)
    {
        return false;
    }

    // Get surrounding faces
    auto t_FacesAroundFace = Helper::init_neighbor_faces(m_Mesh, a_FaceHandle);

    // Check 10 surrounding faces are all quad.
    if(!Helper::are_faces_all_quads(m_Mesh, t_FacesAroundFace))
    {
        return false;
    }

    return true;
}



std::vector<Patch> T2PatchConstructor::getPatch(const FaceHandle& a_FaceHandle)
{
    // Get neighbor verts
    auto t_NBVerts = initNeighborVerts(a_FaceHandle);

    return getPatch(t_NBVerts);
}

std::vector<Patch> T2PatchConstructor::getPatch(const std::vector<VertexHandle>& a_NBVertexHandles)
{
    // Convert Neighbor Verts to matrix type
    auto t_NBVertsMat = Helper::verthandles_to_points_mat(m_Mesh, a_NBVertexHandles);

    // Generate patch
    auto t_BBcoefs = m_Mask * t_NBVertsMat;

    const int a_NumOfPatch = 16;
    return points_mat_to_patches(a_NumOfPatch, "Group 4 T2", t_BBcoefs);
}




std::vector<VertexHandle> T2PatchConstructor::initNeighborVerts(const FaceHandle& a_FaceHandle)
{
    // Init vector for neighbor points
    const int t_NumOfVerts = 20;
    std::vector<VertexHandle> t_NBVerts;
    Helper::set_vert_vector_to_default(t_NumOfVerts, t_NBVerts);

    // Find 9 11 12 13 6 5
    auto t_HexVerts = Helper::get_verts_of_face(m_Mesh, a_FaceHandle);

    // Rotate verts of faces untill the order match t_HexVertOrder
    for(int i=0; i<t_HexVerts.size(); i++)
    {
        bool t_IsCurrentVert3Valence = Helper::is_vert_3_valence(m_Mesh, t_HexVerts[0]);
        bool t_IsNextSecondVert3Valence = Helper::is_vert_3_valence(m_Mesh, t_HexVerts[2]);
        if(t_IsCurrentVert3Valence && t_IsNextSecondVert3Valence)
        {
            break;
        }

        std::rotate(t_HexVerts.begin(), t_HexVerts.begin()+1, t_HexVerts.end());
    }

    // Store 9 11 12 13 6 5
    std::vector<int> t_HexVertOrder = {9,11,12,13,6,5};
    for(int i=0; i<t_HexVerts.size(); i++)
    {
        t_NBVerts[t_HexVertOrder[i]] = t_HexVerts[i];
    }

    // Get the halfedge pointing 6 -> 5
    HalfedgeHandle t_HEHandle;
    for(auto FHIt=m_Mesh.cfh_ccwiter(a_FaceHandle); FHIt.is_valid(); ++FHIt)
    {
        // Get starting point of halfedge (6)
        auto t_StartPointOfHalfedge = m_Mesh.from_vertex_handle(*FHIt);
        if(t_StartPointOfHalfedge==t_NBVerts[6])
        {
            t_HEHandle = *FHIt;
            break;
        }
    }

    // Move halfedge from 6->5 to 5->1
    HalfedgeOperation::get_vert(m_Mesh, t_HEHandle, {3,2,3});

    std::vector<std::vector<int>> t_GetVertOrder = {{1,0},{4},{8},{10,15},{16},{17},{18,19},{14},{7,3},{2}};
    for(int i=0; i<t_GetVertOrder.size(); i++)
    {
        // Pair (2 verts) = corner
        if(t_GetVertOrder[i].size()==2)
        {
            auto t_PairVerts = HalfedgeOperation::get_verts(m_Mesh, t_HEHandle, {1,4,1,4,1,3});
            t_NBVerts[t_GetVertOrder[i][0]] = t_PairVerts[0];
            t_NBVerts[t_GetVertOrder[i][1]] = t_PairVerts[1];
        }
        else
        {
            auto t_Vert = HalfedgeOperation::get_vert(m_Mesh, t_HEHandle, {1,4,1,3});
            t_NBVerts[t_GetVertOrder[i][0]] = t_Vert;
        }
    }



    return t_NBVerts;
}
