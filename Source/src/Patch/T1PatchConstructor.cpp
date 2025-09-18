/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "T1PatchConstructor.hpp"
#include "../Helper/HalfedgeOperation.hpp"
#include "../Helper/ReadCSV2Matrix.hpp"
#include "../Patch/Patch.hpp"

Mat128x18d T1PatchConstructor::getMask()
{
    std::string t_MaskCSVFilePathT1 = "./Table/T1.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathT1, 128, 18); //  NumOfTotalCpts = 128, NumOfVerts = 18
}

bool T1PatchConstructor::isSamePatchType(const FaceHandle& a_FaceHandle)
{
    // Check is the face a Pentagon
    if(!Helper::is_pentagon(m_Mesh, a_FaceHandle))
    {
        return false;
    }

    // Check Pentagon has 4 4D verts and only 1 3D vert
    if(!isPentagonTjunction(a_FaceHandle))
    {
        return false;
    }

    // Check there are 9 neighbor faces around the pentagon
    auto t_FirstLayerFaces = Helper::init_neighbor_faces(m_Mesh, a_FaceHandle);
    if(!Helper::has_9_neighbor_faces(t_FirstLayerFaces))
    {
        return false;
    }

    // Check whether neighbor faces are all quad
    if(!Helper::are_faces_all_quads(m_Mesh, t_FirstLayerFaces))
    {
        return false;
    }

    return true;
}


/*
 * Check if the pentagon only has one tjunction
 */
bool T1PatchConstructor::isPentagonTjunction(const FaceHandle& a_FaceHandle)
{
    const int t_NumOf4ValenceVertsForPenTJ = 4;
    const int t_NumOf3ValenceVertsForPenTJ = 1;
    int t_NumOf4ValenceVerts = 0;
    int t_NumOf3ValenceVerts = 0;
    for (auto t_FVIt = m_Mesh.cfv_iter(a_FaceHandle); t_FVIt.is_valid(); ++t_FVIt)
    {
        if (Helper::is_vert_4_valence(m_Mesh, *t_FVIt))
        {
            t_NumOf4ValenceVerts++;
        }
        else if (Helper::is_vert_3_valence(m_Mesh, *t_FVIt))
        {
            t_NumOf3ValenceVerts++;
        }
    }
    return (t_NumOf4ValenceVertsForPenTJ == t_NumOf4ValenceVerts
            && t_NumOf3ValenceVertsForPenTJ == t_NumOf3ValenceVerts);
}


std::vector<Patch> T1PatchConstructor::getPatch(const FaceHandle& a_FaceHandle)
{
    // Get neighbor verts
    auto t_NBVerts = initNeighborVerts(a_FaceHandle);

    return getPatch(t_NBVerts);
}

std::vector<Patch> T1PatchConstructor::getPatch(const std::vector<VertexHandle>& a_NBVertexHandles)
{
    // Convert Neighbor Verts to matrix type
    auto t_NBVertsMat = Helper::verthandles_to_points_mat(m_Mesh, a_NBVertexHandles);

    // Generate patch
    auto t_BBcoefs = m_Mask * t_NBVertsMat;

    const int a_NumOfPatch = 8;
    return points_mat_to_patches(a_NumOfPatch, "Group 3 T1", t_BBcoefs);
}


/*
 * Get the neighbor verts
 *      0 - 1 - 2 - 3
 *      |   |   |   |
 *      4 - 5 - 6 - 7
 *     /   /     \   \
 *    8 - 9 -10- 11 - 12
 *    |   |   |   |   |
 *   13 -14 -15- 16 - 17
 */
std::vector<VertexHandle> T1PatchConstructor::initNeighborVerts(const FaceHandle& a_FaceHandle)
{
    // Init vector for neighbor points
    const int t_NumOfVerts = 18;
    std::vector<VertexHandle> t_NBVerts;
    Helper::set_vert_vector_to_default(t_NumOfVerts, t_NBVerts);

    // Find 10 11 6 5 9  **Verts from getVertsOfFace are in clockwise order
    auto t_PenVerts = Helper::get_verts_of_face(m_Mesh, a_FaceHandle);
    std::vector<int> t_PenVertOrder{10, 9, 5, 6, 11};
    while(true)
    {
        //rotate penverts until its first element is 3D vert
        if(Helper::is_vert_3_valence(m_Mesh, t_PenVerts[0]))
        {
            break;
        }
        std::rotate(t_PenVerts.begin(),t_PenVerts.begin()+1,t_PenVerts.end());
    }
    for(int i=0; i<t_PenVerts.size(); i++)
    {
        t_NBVerts[t_PenVertOrder[i]] = t_PenVerts[i];
    }

    // Find 15 and the halfedge (10->15)
    MeshType::HalfedgeHandle t_CurrentHef;
    for(auto t_VHIt = m_Mesh.cvoh_iter(t_NBVerts[10]); t_VHIt.is_valid(); ++t_VHIt)
    {
        auto t_NextHalfedge = m_Mesh.next_halfedge_handle(*t_VHIt);
        auto t_Vert = m_Mesh.from_vertex_handle(t_NextHalfedge);
        bool t_IsVertNotInPen = std::find(t_PenVerts.begin(), t_PenVerts.end(), t_Vert) != t_PenVerts.end();
        if(!t_IsVertNotInPen)
        {
            t_NBVerts[15] = t_Vert;
            t_CurrentHef = *t_VHIt;
            break;
        }
    }

    // Find 4 corners 14 13 18 -> 4 0 1 -> 2 3 7 -> 12 17 16
    std::vector<std::vector<int>> t_GetVertOrder{{14, 13, 8}, {4, 0, 1}, {2, 3, 7}, {12, 17, 16}};
    // Operation for one corner. Operations: 1=next 2=previous 3=opposite 4=get
    std::vector<int> t_CornerOperation{1, 1, 4, 3, 1, 1, 4, 1, 4, 3};

    HalfedgeOperation::init_verts(m_Mesh, t_CurrentHef, t_NBVerts, t_GetVertOrder, t_CornerOperation);

    return t_NBVerts;
}
