/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "T0PatchConstructor.hpp"
#include "../Helper/HalfedgeOperation.hpp"
#include "../Helper/ReadCSV2Matrix.hpp"
#include "../Patch/Patch.hpp"

Mat64x14d T0PatchConstructor::getMask()
{
    std::string t_MaskCSVFilePathT0 = "./Table/T0.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathT0, 64, 14); //  NumOfTotalCpts = 64, NumOfVerts = 14
}

bool T0PatchConstructor::isSamePatchType(const FaceHandle& a_FaceHandle)
{
    // The face should be triangle
    if(!Helper::is_triangle(m_Mesh, a_FaceHandle))
    {
        return false;
    }

    // The verts of triangle should be one 5 valences and two 4 valences
    int t_5valenceVertCounter = 0;
    int t_4valenceVertCounter = 0;
    for(auto FVIt = m_Mesh.cfv_ccwiter(a_FaceHandle); FVIt.is_valid(); ++FVIt)
    {
        if(Helper::is_vert_5_valence(m_Mesh, *FVIt))
        {
            t_5valenceVertCounter++;
        }
        else if(Helper::is_vert_4_valence(m_Mesh, *FVIt))
        {
            t_4valenceVertCounter++;
        }
    }
    if(t_5valenceVertCounter!=1 || t_4valenceVertCounter!=2)
    {
        return false;
    }

    // There should be 7 faces around the triangle as first layer
    auto t_FirstLayerFaces = Helper::init_neighbor_faces(m_Mesh, a_FaceHandle);
    if(!Helper::has_7_neighbor_faces(t_FirstLayerFaces))
    {
        return false;
    }

    // First layer faces should be all quads
    if(!Helper::are_faces_all_quads(m_Mesh, t_FirstLayerFaces))
    {
        return false;
    }
    
    return true;
}

std::vector<Patch> T0PatchConstructor::getPatch(const FaceHandle& a_FaceHandle)
{
    // Get neighbor verts
    auto t_NBVerts = initNeighborVerts(a_FaceHandle);

    return getPatch(t_NBVerts);
}

std::vector<Patch> T0PatchConstructor::getPatch(const std::vector<VertexHandle>& a_NBVertexHandles)
{
    // Convert Neighbor Verts to matrix type
    auto t_NBVertsMat = Helper::verthandles_to_points_mat(m_Mesh, a_NBVertexHandles);

    // Generate patch
    auto t_BBcoefs = m_Mask * t_NBVertsMat;

    const int a_NumOfPatch = 4;
    return points_mat_to_patches(a_NumOfPatch, "Group 2 T0", t_BBcoefs);
}


/*
 * Get the neighbor verts
 *      0 - 1 - 2
 *      |   |   |
 *      3 - 4 - 5
 *     /   / \   \
 *    6 - 7 - 8 - 9
 *    |   |   |   |
 *   10 -11 -12 - 13
 */
std::vector<VertexHandle> T0PatchConstructor::initNeighborVerts(const FaceHandle& a_FaceHandle)
{
    // Init vector for neighbor points
    const int t_NumOfVerts = 14;
    std::vector<VertexHandle> t_NBVerts;
    Helper::set_vert_vector_to_default(t_NumOfVerts, t_NBVerts);

    // Get triangle verts 4 7 8
    std::vector<VertexHandle> t_TriangleVerts;
    for(auto FVIt = m_Mesh.cfv_ccwiter(a_FaceHandle); FVIt.is_valid(); ++FVIt)
    {
        t_TriangleVerts.push_back(*FVIt);
    }

    //Rotate vector to make five valence vert at postion 0
    while(true)
    {
        if(Helper::is_vert_5_valence(m_Mesh, t_TriangleVerts[0]))
        {
            break;
        }
        std::rotate(t_TriangleVerts.begin(), t_TriangleVerts.begin()+1, t_TriangleVerts.end());
    }

    // Store the verts into array
    std::vector<int> t_TriVertOrder = {4, 7, 8};
    for(int i=0; i<t_TriVertOrder.size(); i++)
    {
        t_NBVerts[t_TriVertOrder[i]] = t_TriangleVerts[i];
    }

    // Get the halfedge 8->4
    HalfedgeHandle t_HalfedgeHandle;
    for(auto FHIt=m_Mesh.cfh_ccwiter(a_FaceHandle); FHIt.is_valid(); ++FHIt)
    {
        if(t_NBVerts[8] == m_Mesh.from_vertex_handle(*FHIt))
        {
            t_HalfedgeHandle = *FHIt;
            break;
        }
    }

    // Get vert 5
    auto t_VertNo5 = HalfedgeOperation::get_verts(m_Mesh, t_HalfedgeHandle, {3, 2, 4});
    t_NBVerts[5] = t_VertNo5[0];

    // Get verts 2 1 0 3
    auto t_TopCommands = HalfedgeOperation::duplicate_commands(2, {3, 1, 1, 4, 1, 4});
    std::vector<int> t_TopVertOrder = {2, 1, 0, 3};
    HalfedgeOperation::init_verts(m_Mesh, t_HalfedgeHandle, t_NBVerts, t_TopVertOrder, t_TopCommands);

    // Get vert 6 10 11 12 13 9
    auto t_BottomCommands = HalfedgeOperation::duplicate_commands(2, {3, 1, 1, 4, 3, 1, 1, 4, 1, 4});
    std::vector<int> t_BottomVertOrder = {6, 10, 11, 12, 13, 9};
    HalfedgeOperation::init_verts(m_Mesh, t_HalfedgeHandle, t_NBVerts, t_BottomVertOrder, t_BottomCommands);

    return t_NBVerts;
}
