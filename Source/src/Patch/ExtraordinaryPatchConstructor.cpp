/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "ExtraordinaryPatchConstructor.hpp"
#include "../Helper/HalfedgeOperation.hpp"
#include "../Helper/ReadCSV2Matrix.hpp"
#include "../Patch/Patch.hpp"

Matrix ExtraordinaryPatchConstructor::getMask()
{
    switch(m_Valence)
    {
        case 3:
            return m_MaskSct3;
        case 5:
            return m_MaskSct5;
        case 6:
            return m_MaskSct6;
        case 7:
            return m_MaskSct7;
        case 8:
            return m_MaskSct8;
    }
    std::cerr << "ERROR: Extraordinary point valence must be either 3, 5, 6, 7, or 8" << std::endl;
    return Matrix();
}

Mat48x7d ExtraordinaryPatchConstructor::getMaskSct3()
{
    std::string t_MaskCSVFilePathSct3 = "./Table/eopSct3.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct3, 48, 7);
}

Mat80x11d ExtraordinaryPatchConstructor::getMaskSct5()
{
    std::string t_MaskCSVFilePathSct5 = "./Table/eopSct5.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct5, 80, 11);
}

Mat384x13d ExtraordinaryPatchConstructor::getMaskSct6()
{
    std::string t_MaskCSVFilePathSct6 = "./Table/eopSct6.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct6, 384, 13);
}

Mat448x15d ExtraordinaryPatchConstructor::getMaskSct7()
{
    std::string t_MaskCSVFilePathSct7 = "./Table/eopSct7.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct7, 448, 15);
}

Mat512x17d ExtraordinaryPatchConstructor::getMaskSct8()
{
    std::string t_MaskCSVFilePathSct8 = "./Table/eopSct8.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct8, 512, 17);
}

bool ExtraordinaryPatchConstructor::isSamePatchType(const VertexHandle& a_VertexHandle)
{
    // Extraordinary point should be 3, 5, 6, 7, or 8 valence
    m_Valence = Helper::get_vert_valence(m_Mesh, a_VertexHandle);
    if(m_Valence<3 || m_Valence==4 || m_Valence>8)
    {
        return false;
    }

    // Can not be on the edge
    if(m_Mesh.is_boundary(a_VertexHandle))
    {
        return false;
    }

    // // The first layer of faces should all be quad
    // auto t_FirstLayerFaces = Helper::get_faces_around_vert_counterclock(m_Mesh, a_VertexHandle);
    // for(auto t_Vert : t_FirstLayerFaces)
    // {
    //     if(!Helper::are_faces_all_quads(m_Mesh, t_FirstLayerFaces))
    //     {
    //         return false;
    //     }
    // }

    // The surrounding faces should be either quad or triangle and no more than two
    // consecutive triangles.
    auto t_FHs = Helper::get_faces_around_vert_counterclock(m_Mesh, a_VertexHandle);
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
 * The index of neighbor verts (following example is for 3 valent point)
 *        0
 *      / | \
 *     /  |  \
 *    1   6   5
 *    \  / \ /
 *     2    4
 *      \  /
 *       3
 */
std::vector<VertexHandle> ExtraordinaryPatchConstructor::initNeighborVerts(const VertexHandle& a_VertexHandle)
{
    std::vector<VertexHandle> t_VertexHandles;
    std::vector<VertexHandle> t_NBVerts;
    for(auto t_VHIt = m_Mesh.cvoh_ccwiter(a_VertexHandle); t_VHIt.is_valid(); ++t_VHIt)
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
            t_VertexHandles.push_back(t_NBVert);
        }
    }

    // store the central point
    t_VertexHandles.push_back(a_VertexHandle);

    return t_VertexHandles;
}


std::vector<Patch> ExtraordinaryPatchConstructor::getPatch(const VertexHandle& a_VertexHandle)
{
    // Get neighbor verts
    auto t_NBVerts = initNeighborVerts(a_VertexHandle);

    return getPatch(t_NBVerts);
}

std::vector<Patch> ExtraordinaryPatchConstructor::getPatch(const std::vector<VertexHandle>& a_NBVertexHandles)
{
    // Convert neighbor verts to matrix type
    auto t_NBVertsMat = Helper::verthandles_to_points_mat(m_Mesh, a_NBVertexHandles);

    // multiply CC points with mask to generate patches
    Matrix t_BBcoefs = getMask() * t_NBVertsMat;

    int t_NumOfPatches = m_Valence;
    if(m_Valence >= 6)
    {
        t_NumOfPatches *= 4; // have four patches in each sector
    }


    // Save control points into Patch
    return points_mat_to_patches(t_NumOfPatches, "Group 1 ExtraordinaryPoint", t_BBcoefs);
}
