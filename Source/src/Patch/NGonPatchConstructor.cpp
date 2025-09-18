/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "NGonPatchConstructor.hpp"
#include "../Helper/HalfedgeOperation.hpp"
#include "../Helper/ReadCSV2Matrix.hpp"
#include "../Patch/Patch.hpp"

Matrix NGonPatchConstructor::getMask()
{
    switch(m_FaceValence)
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
    std::cerr << "ERROR: N-Gon face must have face-valence either 3, 5, 6, 7, or 8" << std::endl;
    return Matrix();
}

Mat48x12d NGonPatchConstructor::getMaskSct3()
{
    std::string t_MaskCSVFilePathSct3 = "./Table/ngonSct3.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct3, 48, 12);
}

Mat80x20d NGonPatchConstructor::getMaskSct5()
{
    std::string t_MaskCSVFilePathSct5 = "./Table/ngonSct5.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct5, 80, 20);
}

Mat384x24d NGonPatchConstructor::getMaskSct6()
{
    std::string t_MaskCSVFilePathSct6 = "./Table/ngonSct6.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct6, 384, 24);
}

Mat448x28d NGonPatchConstructor::getMaskSct7()
{
    std::string t_MaskCSVFilePathSct7 = "./Table/ngonSct7.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct7, 448, 28);
}

Mat512x32d NGonPatchConstructor::getMaskSct8()
{
    std::string t_MaskCSVFilePathSct8 = "./Table/ngonSct8.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct8, 512, 32);
}

bool NGonPatchConstructor::isSamePatchType(const FaceHandle& a_FaceHandle)
{
    // The face should have 3 5 6 7 8 vertices
    m_FaceValence = Helper::get_num_of_verts_for_face(m_Mesh, a_FaceHandle);
    if(m_FaceValence < 3 || m_FaceValence > 8 || m_FaceValence == 4)
    {
        return false;
    }

    // All the vertices of central face should be 4-valent
    if(!Helper::are_verts_of_face_all_4_valence(m_Mesh, a_FaceHandle))
    {
        return false;
    }

    // All surrounded faces should be quads
    auto t_FirstLayerFaces = Helper::init_neighbor_faces(m_Mesh, a_FaceHandle);
    if(!Helper::are_faces_all_quads(m_Mesh, t_FirstLayerFaces))
    {
        return false;
    }

    return true;
}

std::vector<Patch> NGonPatchConstructor::getPatch(const FaceHandle& a_FaceHandle)
{
    // Get neighbor verts
    auto t_NBVerts = initNeighborVerts(a_FaceHandle);

    return getPatch(t_NBVerts);
}

std::vector<Patch> NGonPatchConstructor::getPatch(const std::vector<VertexHandle>& a_NBVertexHandles)
{
    // Convert Neighbor Verts to matrix type
    auto t_NBVertsMat = Helper::verthandles_to_points_mat(m_Mesh, a_NBVertexHandles);

    // Generate patch
    Matrix t_BBcoefs = getMask() * t_NBVertsMat;

    // Only one bi3 patch for each sector when n = 3,5 other wise four per sector
    int a_NumOfPatch = m_FaceValence;
    if(m_FaceValence!=3 && m_FaceValence!=5)
    {
        a_NumOfPatch = a_NumOfPatch * 4;
    }

    return points_mat_to_patches(a_NumOfPatch, "Group 8 nGon", t_BBcoefs);
}


/*
 * Get the neighbor verts
 *           8
 *        /    \
 *      9 - 11- 10
 *     /   / \   \
 *    2 - 3 - 7 - 5
 *    |   |   |   |
 *    0 - 1 - 6 - 4
 */
std::vector<VertexHandle> NGonPatchConstructor::initNeighborVerts(const FaceHandle& a_FaceHandle)
{
    // Init vector for neighbor points
    const int t_NumOfVerts = 4 * m_FaceValence;
    std::vector<VertexHandle> t_NBVerts;
    Helper::set_vert_vector_to_default(t_NumOfVerts, t_NBVerts);

    // get vertices corner by corner (0 1 2 3 is a corner)
    // Get the halfedge 3->7
    int t_WingID = 0;
    for(auto t_FHIt = m_Mesh.cfh_ccwiter(a_FaceHandle); t_FHIt.is_valid(); t_FHIt++)
    {
        // Move halfedge to 0->1 and get vert 0 1 3 2
        int k = t_WingID * 4;
        std::vector<int> t_IDs = {3+k, 1+k, 0+k, 2+k};
        HalfedgeOperation::init_verts_fixed_halfedge(m_Mesh, *t_FHIt, t_NBVerts, t_IDs, {3,1,4,3,4,2,4,2,4});

        t_WingID++;
    }

    return t_NBVerts;
}
