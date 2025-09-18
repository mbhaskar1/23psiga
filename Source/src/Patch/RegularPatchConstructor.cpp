/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "RegularPatchConstructor.hpp"
#include "../Helper/HalfedgeOperation.hpp"
#include "../Patch/Patch.hpp"

Matrix RegularPatchConstructor::getMask(){
    return m_Mask;
}

/*
 * Check if the current face (facehandle) and its neighbors match the Regular structure
 */
bool RegularPatchConstructor::isSamePatchType(const VertexHandle& a_VertHandle)
{
    // Check if the vertex four valence
    if(!Helper::is_vert_4_valence(m_Mesh, a_VertHandle))
    {
        return false;
    }

    // The surrounding faces should be either quad or triangle and no more than two
    // consecutive triangles.
    auto t_FHs = Helper::get_faces_around_vert_counterclock(m_Mesh, a_VertHandle);
    bool t_IsPrevTri = false;
    for(int i=0; i<t_FHs.size(); i++)
    {
        if(!Helper::is_triangle(m_Mesh, t_FHs[i]) && !Helper::is_quad(m_Mesh, t_FHs[i]))
        {
            return false;
        }

        if(Helper::is_triangle(m_Mesh, t_FHs[i]) && 
           Helper::is_triangle(m_Mesh, t_FHs[(i+1)%t_FHs.size()]) &&
           Helper::is_triangle(m_Mesh, t_FHs[(i+2)%t_FHs.size()]))
        {
            return false;
        }
    }

    return true;
}


/*
 * Find the neighbor verts around the given face, then use the verts to generate patch
 */
std::vector<Patch> RegularPatchConstructor::getPatch(const VertexHandle& a_VertHandle)
{
    auto t_NBVertexHandles = initNeighborVerts(a_VertHandle);
    return getPatch(t_NBVertexHandles);
}

/*
 * Given neighbor verts, multiply it with mask to generate patch.
 */
std::vector<Patch> RegularPatchConstructor::getPatch(const std::vector<VertexHandle>& a_NBVertexHandles)
{
    // Convert Neighbor Verts to matrix type
    auto t_NBVertsMat = Helper::verthandles_to_points_mat(m_Mesh, a_NBVertexHandles);

    // Generate patch
    auto t_BBcoefs = m_Mask * t_NBVertsMat;

    const int t_PatchDegU = 2;
    const int t_PatchDegV = 2;
    std::vector<Patch> t_Patches;
    t_Patches.push_back(points_mat_to_patch(t_PatchDegU, t_PatchDegV, "Group 0 Regular", t_BBcoefs));
    return t_Patches;
}


/*
 * Initialize (obtain) the vertices on original mesh needed for getPatch
 */
std::vector<VertexHandle> RegularPatchConstructor::initNeighborVerts(const VertexHandle& a_VertHandle)
{
    // Init m_NeighborVerts with 9 default vertexhandle
    const int t_NumOfVerts = 9;
    std::vector<VertexHandle> t_NBVertexHandles;
    Helper::set_vert_vector_to_default(t_NumOfVerts, t_NBVertexHandles);

    // Get center vertex
    t_NBVertexHandles[4] = a_VertHandle;

    // Get the halfedge belonds to current face
    auto t_CurrentHef = m_Mesh.halfedge_handle(a_VertHandle);

    /* Using halfedge operation to get value for each neighbors
     * Operations: 1=next 2=previous 3=opposite 4=get
     *  Vertex index:
     *    0  1  2
     *    3  4  5
     *    6  7  8
     */
    // Four wings
    std::vector<int> t_GetVertOrder = {1,0,3,6,7,8,5,2};

    // // Operation for one wing
    // std::vector<int> t_WingOperation{1,4,1,4,1,3};

    // HalfedgeOperation::init_verts(m_Mesh, t_CurrentHef, t_NBVertexHandles, t_GetVertOrder, t_WingOperation);
    int i = 0;
    std::vector<VertexHandle> t_NBVerts;
    for(auto t_VHIt = m_Mesh.cvoh_ccwiter(a_VertHandle); t_VHIt.is_valid(); ++t_VHIt)
    {
        auto t_CurrFH = m_Mesh.face_handle(*t_VHIt);
        std::vector<int> t_Commands;
        if(Helper::is_quad(m_Mesh, t_CurrFH))
        {
            t_Commands = {1,4,1,4};
        }
        else if(Helper::is_triangle(m_Mesh, t_CurrFH))
        {   
            auto t_OppositeHE = m_Mesh.opposite_halfedge_handle(*t_VHIt);
            auto t_PreFH = m_Mesh.face_handle(t_OppositeHE);
            if(Helper::is_quad(m_Mesh, t_PreFH))
            {
                t_Commands = {1,4,1,4};
            }
            else
            {
                t_Commands = {1,4,4};
            }
        }
        else
        {
            assert(false);
        }

        t_NBVerts = HalfedgeOperation::get_verts_fixed_halfedge(m_Mesh, *t_VHIt, t_Commands);
        
        for(auto t_NBVert : t_NBVerts)
        {
            t_NBVertexHandles[t_GetVertOrder[i]] = t_NBVert;
            i++;
        }
    }

    return t_NBVertexHandles;
}
