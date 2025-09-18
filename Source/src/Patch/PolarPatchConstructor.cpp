/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "PolarPatchConstructor.hpp"
#include "../Helper/HalfedgeOperation.hpp"
#include "../Helper/ReadCSV2Matrix.hpp"
#include "../Patch/Patch.hpp"

Matrix PolarPatchConstructor::getMask(){
    switch (m_NumOfSct)
    {
        case 3:
            return m_MaskSct3;
        case 4:
            return m_MaskSct4;
        case 5:
            return m_MaskSct5;
        case 6:
            return m_MaskSct6;
        case 7:
            return m_MaskSct7;
        case 8:
            return m_MaskSct8;
    }
    std::cerr << "ERROR: Polar point must have 3 to 8 valence" << std::endl;
    return Matrix();
}

/*
 * Get the mask for generating Bi3 Patch
 */
Mat36x4d PolarPatchConstructor::getMaskSct3()
{
    std::string t_MaskCSVFilePathSct3 = "./Table/polarSct3.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct3, 36, 4); //  NumOfTotalCpts = 60, NumOfVerts = 4
}

Mat48x5d PolarPatchConstructor::getMaskSct4()
{
    std::string t_MaskCSVFilePathSct4 = "./Table/polarSct4.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct4, 48, 5);
}

Mat60x6d PolarPatchConstructor::getMaskSct5()
{
    std::string t_MaskCSVFilePathSct5 = "./Table/polarSct5.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct5, 60, 6);
}

Mat72x7d PolarPatchConstructor::getMaskSct6()
{
    std::string t_MaskCSVFilePathSct6 = "./Table/polarSct6.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct6, 72, 7);
}

Mat84x8d PolarPatchConstructor::getMaskSct7()
{
    std::string t_MaskCSVFilePathSct7 = "./Table/polarSct7.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct7, 84, 8);
}

Mat96x9d PolarPatchConstructor::getMaskSct8()
{
    std::string t_MaskCSVFilePathSct8 = "./Table/polarSct8.csv";
    return read_csv_as_matrix(t_MaskCSVFilePathSct8, 96, 9);
}


bool PolarPatchConstructor::isSamePatchType(const VertexHandle& a_VertexHandle)
{
    // Polar point should be 3 to 8 valence
    m_NumOfSct = Helper::get_vert_valence(m_Mesh, a_VertexHandle);
    const int t_Lowerbound = 3;
    const int t_Upperbound = 8;
    if(m_NumOfSct < t_Lowerbound || m_NumOfSct > t_Upperbound)
    {
        return false;
    }

    // Polar point should not be on the boundary
    if(m_Mesh.is_boundary(a_VertexHandle))
    {
        return false;
    }

    // The first layer of surrounded faces should all be triangles
    auto t_NBFaces = Helper::get_faces_around_vert_counterclock(m_Mesh, a_VertexHandle);
    for(auto t_Face : t_NBFaces)
    {
        if(!Helper::is_triangle(m_Mesh, t_Face))
        {
            return false;
        }
    }

    return true;
}


std::vector<Patch> PolarPatchConstructor::getPatch(const VertexHandle& a_VertexHandle)
{
    auto t_NBVerts = initNeighborVerts(a_VertexHandle);

    return getPatch(t_NBVerts);
}

std::vector<Patch> PolarPatchConstructor::getPatch(const std::vector<VertexHandle>& a_NBVertexHandles)
{
    // Convert NeighborVerts to matrix type
    auto t_NBVertsMat = Helper::verthandles_to_points_mat(m_Mesh, a_NBVertexHandles);

    // Get mask
    Matrix t_BBcoefs = getMask() * t_NBVertsMat;

    const int t_DegU = 3;
    const int t_DegV = 2;
    auto t_Patches = points_mat_to_patches(t_DegU, t_DegV, "Group 5 Polar", t_BBcoefs);

    return t_Patches;
}


/*
 *     P1 ---- P4
 *      | \  / |
 *      |  P0  |
 *      | /  \ |
 *     P2 ---- P3
 */
std::vector<VertexHandle> PolarPatchConstructor::initNeighborVerts(const VertexHandle& a_VertexHandle)
{
    // init vector to stroe all vertices
    std::vector<VertexHandle> t_NBVertexHandles;

    // Get first layer vetices  ex: for 4sct, P1 -> P2 -> P3 -> P4
    for(auto VHIt=m_Mesh.cvoh_ccwiter(a_VertexHandle); VHIt.is_valid(); ++VHIt)
    {
        t_NBVertexHandles.push_back(HalfedgeOperation::get_vert_fixed_halfedge(m_Mesh, *VHIt, {1,4}));
    }

    // Insert the central point P0
    t_NBVertexHandles.insert(t_NBVertexHandles.begin(), a_VertexHandle); ;

    return t_NBVertexHandles;
}
