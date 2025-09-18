/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "Helper.hpp"
#include "HalfedgeOperation.hpp"
#include <iostream>
#include <cmath>
#include <unordered_set>

namespace Helper
{

// Vert functions

/**
 * Get the valence of a vertex in a mesh, i.e., the number of edges connected to the vertex
 * @param a_Mesh the given mesh
 * @param a_VertexHandle the given vertex
 * @return the valence of the given vertex in the given mesh
 */
int get_vert_valence(const MeshType& a_Mesh, const VertexHandle& a_VertexHandle)
{
    return a_Mesh.valence(a_VertexHandle);
}

/**
 * Checks if a vertex in a mesh has valence 3 (i.e., if the vertex has 3 edges connected to it)
 * @param a_Mesh the given mesh
 * @param a_VertexHandle the given vertex
 * @return true if the valence of the vertex is 3; false otherwise
 */
bool is_vert_3_valence(const MeshType& a_Mesh, const VertexHandle& a_VertexHandle)
{
    return (get_vert_valence(a_Mesh, a_VertexHandle)==3) ? true : false;
}

/**
 * Checks if a vertex in a mesh has valence 4 (i.e., if the vertex has 4 edges connected to it)
 * @param a_Mesh the given mesh
 * @param a_VertexHandle the given vertex
 * @return true if the valence of the vertex is 4; false otherwise
 */
bool is_vert_4_valence(const MeshType& a_Mesh, const VertexHandle& a_VertexHandle)
{
    return (get_vert_valence(a_Mesh, a_VertexHandle)==4) ? true : false;
}

/**
 * Checks if a vertex in a mesh has valence 5 (i.e., if the vertex has 5 edges connected to it)
 * @param a_Mesh the given mesh
 * @param a_VertexHandle the given vertex
 * @return true if the valence of the vertex is 5; false otherwise
 */
bool is_vert_5_valence(const MeshType& a_Mesh, const VertexHandle& a_VertexHandle)
{
    return (get_vert_valence(a_Mesh, a_VertexHandle)==5) ? true : false;
}

/**
 * Checks if all vertices of a face in a mesh have valence 4 (i.e., if each vertex has 4 edges connected to it).
 * @param a_Mesh the given mesh
 * @param a_FaceHandle the given face
 * @return true if the valence of each vertex of the given face is 4; false otherwise
 */
bool are_verts_of_face_all_4_valence(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle)
{
    for(auto FVIt = a_Mesh.cfv_iter(a_FaceHandle); FVIt.is_valid(); FVIt++)
    {
        if(!is_vert_4_valence(a_Mesh, *FVIt))
        {
            return false;
        }
    }
    return true;
}

void set_vert_vector_to_default(const int a_Size, std::vector<VertexHandle>& a_VertexHandles)
{
    a_VertexHandles.clear();
    VertexHandle t_DefaultVertexHandle;
    for(int i=0; i<a_Size; i++)
    {
        a_VertexHandles.push_back(t_DefaultVertexHandle);
    }
}

// Return all the facehandles around the given vertexhandle
std::vector<FaceHandle> get_faces_around_vert_counterclock(const MeshType& a_Mesh, const VertexHandle& a_VertHandle)
{
    std::vector<FaceHandle> t_Faces;
    for(auto t_VFIt = a_Mesh.cvf_ccwiter(a_VertHandle); t_VFIt.is_valid(); ++t_VFIt)
    {
        t_Faces.push_back(*t_VFIt);
    }
    return t_Faces;
}

/*
 *  Get first and second layers of faces aroumnd the vert (unorder)
 *  ex:
 *      o - o - o - o - o
 *      |   |   |   |   |
 *      o - o - o - o - o
 *      |   |   |   |   |
 *      o - o - v - o - o
 *      |   |   |   |   |
 *      o - o - o - o - o
 *      |   |   |   |   |
 *      o - o - o - o - o
 *
 */
std::vector<FaceHandle> get_two_layers_faces_around_vert(const MeshType& a_Mesh, const VertexHandle& a_VertHandle)
{
    std::vector<FaceHandle> t_TwoLayerFaces;

    // Get first layer of faces around vertex
    // And get the neighbor faces of each first layer face (neighbors might be duplicated)
    for(auto VF_It = a_Mesh.cvf_iter(a_VertHandle); VF_It.is_valid(); ++VF_It)
    {
        auto t_FaceNBFaces = init_neighbor_faces(a_Mesh, *VF_It);
        t_TwoLayerFaces.insert(t_TwoLayerFaces.end(), t_FaceNBFaces.begin(), t_FaceNBFaces.end());
    }

    // Eliminate the duplicated faces
    sort(t_TwoLayerFaces.begin(), t_TwoLayerFaces.end());
    t_TwoLayerFaces.erase(unique(t_TwoLayerFaces.begin(), t_TwoLayerFaces.end()), t_TwoLayerFaces.end());

    return t_TwoLayerFaces;
}

/*
 * Return two layer of vertexhandle around a given vert
 * w/o ordering
 */
std::vector<VertexHandle> get_two_layers_verts_around_vert(const MeshType& a_Mesh, const VertexHandle& a_VertHandle)
{
    std::vector<VertexHandle> t_AllVerts;
    auto t_TwoLayersFaces = get_two_layers_faces_around_vert(a_Mesh, a_VertHandle);

    //Get vertices from each face
    for(auto t_Face : t_TwoLayersFaces)
    {
        auto t_VertsForAFace = get_verts_of_face(a_Mesh, t_Face);
        t_AllVerts.insert(t_AllVerts.end(), t_VertsForAFace.begin(), t_VertsForAFace.end());
    }

    // Eliminate the duplicated verts
    sort(t_AllVerts.begin(), t_AllVerts.end());
    t_AllVerts.erase(unique(t_AllVerts.begin(), t_AllVerts.end()), t_AllVerts.end());

    return t_AllVerts;
}




// Face functions
/**
 * Checks if a face is a triangle
 * @param a_Mesh the given mesh
 * @param a_FaceHandle the given face
 * @return true if the face is a triangle; false otherwise
 */
bool is_triangle(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle)
{
    const int t_NumOfVertsForQuad = 3;
    return (get_num_of_verts_for_face(a_Mesh, a_FaceHandle)==t_NumOfVertsForQuad) ? true : false;
}

/**
 * Checks if a face is a quadrilateral
 * @param a_Mesh the given mesh
 * @param a_FaceHandle the given face
 * @return true if the face is a quadrilateral; false otherwise
 */
bool is_quad(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle)
{
    const int t_NumOfVertsForQuad = 4;
    return (get_num_of_verts_for_face(a_Mesh, a_FaceHandle)==t_NumOfVertsForQuad) ? true : false;
}

/**
 * Checks if a face is a pentagon
 * @param a_Mesh the given mesh
 * @param a_FaceHandle the given face
 * @return true if the face is a pentagon; false otherwise
 */
bool is_pentagon(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle)
{
    const int t_NumOfVertsForPentagon = 5;
    return (get_num_of_verts_for_face(a_Mesh, a_FaceHandle)==t_NumOfVertsForPentagon) ? true : false;
}

/**
 * Checks if a face is a hexagon
 * @param a_Mesh the given mesh
 * @param a_FaceHandle the given face
 * @return true if the face is a hexagon; false otherwise
 */
bool is_hexagon(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle)
{
    const int t_NumOfVertsForHexagon = 6;
    return (get_num_of_verts_for_face(a_Mesh, a_FaceHandle)==t_NumOfVertsForHexagon) ? true : false;
}

/**
 * Gets a list containing all neighboring faces of some given face
 * @param a_Mesh the given mesh
 * @param a_FaceHandle the given face
 * @return A vector of FaceHandle's containing handles of each face neighboring the given face in the given mesh
 */
std::vector<FaceHandle> init_neighbor_faces(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle)
{
    std::vector<FaceHandle> t_NBFaceHandles;

    for (auto t_FVIt = a_Mesh.cfv_iter(a_FaceHandle); t_FVIt.is_valid(); ++t_FVIt)
    {
        for (auto t_VFIt = a_Mesh.cvf_iter(*t_FVIt); t_VFIt.is_valid(); ++t_VFIt)
        {
            t_NBFaceHandles.push_back(*t_VFIt);
        }
    }
    // Remove center face (Original)
    t_NBFaceHandles.erase(std::remove(t_NBFaceHandles.begin(),
                          t_NBFaceHandles.end(), a_FaceHandle), t_NBFaceHandles.end());

    // Remove duplicated neighbor faces
    sort(t_NBFaceHandles.begin(), t_NBFaceHandles.end());
         t_NBFaceHandles.erase(unique(t_NBFaceHandles.begin(),
                               t_NBFaceHandles.end()), t_NBFaceHandles.end());

    return t_NBFaceHandles;
}

bool has_7_neighbor_faces(const std::vector<FaceHandle>& a_NBFaceHandles)
{
    const int t_7Neighbor = 7;
    return (get_num_of_neighbor_faces(a_NBFaceHandles)==t_7Neighbor) ? true : false;
}

bool has_9_neighbor_faces(const std::vector<FaceHandle>& a_NBFaceHandles)
{
    const int t_9Neighbor = 9;
    return (get_num_of_neighbor_faces(a_NBFaceHandles)==t_9Neighbor) ? true : false;
}

/**
 * Checks if each face in a list of faces is a quadrilateral
 * @param a_Mesh the given mesh
 * @param a_FaceHandles a vector of faces
 * @return true if each face in the given vector of faces is a quadrilateral; false otherwise
 */
bool are_faces_all_quads(const MeshType& a_Mesh, const std::vector<FaceHandle>& a_FaceHandles)
{
    for(auto t_FH : a_FaceHandles)
    {
        if(!is_quad(a_Mesh, t_FH))
        {
            return false;
        }
    }
    return true;
}

int get_num_of_verts_for_face(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle)
{
    return a_Mesh.valence(a_FaceHandle);
}

int get_num_of_neighbor_faces(const std::vector<FaceHandle>& a_NBFaceHandles)
{
    return a_NBFaceHandles.size();
}

std::vector<VertexHandle> get_verts_of_face(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle)
{
    std::vector<VertexHandle> t_VertsOfFace;
    for (auto t_FVIt = a_Mesh.cfv_ccwiter(a_FaceHandle); t_FVIt.is_valid(); ++t_FVIt)
    {
        t_VertsOfFace.push_back(*t_FVIt);
    }
    return t_VertsOfFace;
}

int num_of_quads(const MeshType& a_Mesh, std::vector<FaceHandle> a_FaceHandles)
{
    int t_NumOfQuads = 0;
    for(auto t_Face : a_FaceHandles)
    {
        if(is_quad(a_Mesh, t_Face))
        {
            t_NumOfQuads++;
        }
    }
    return t_NumOfQuads;
}

int num_of_triangles(const MeshType& a_Mesh, std::vector<FaceHandle> a_FaceHandles)
{
    int t_NumOfTriangles = 0;
    for(auto t_Face : a_FaceHandles)
    {
        if(is_triangle(a_Mesh, t_Face))
        {
            t_NumOfTriangles++;
        }
    }
    return t_NumOfTriangles;
}




// Type conversion

Vec3d verthandles_to_point_vec(const MeshType& a_Mesh, const VertexHandle& a_VertHandle)
{
    Vec3d t_PointVec = {0, 0, 0};
    if(a_VertHandle.is_valid())
    {
        for(int i=0; i<3; i++)
        {
            t_PointVec[i] = a_Mesh.point(a_VertHandle)[i];
        }
    }
    return t_PointVec;
}

Matrix verthandles_to_points_mat(const MeshType& a_Mesh, const std::vector<VertexHandle>& a_VertHandle)
{
    int t_NumOfVerts = a_VertHandle.size();
    auto t_PointMat = Matrix(t_NumOfVerts, 3);

    for(int i=0; i<t_NumOfVerts; i++)
    {
        Vec3d point = verthandles_to_point_vec(a_Mesh, a_VertHandle[i]);
        for(int j=0; j < 3; j++){
            t_PointMat(i, j) = point[j];
        }
    }
    return t_PointMat;
}

// Gauss Point Evaluation
/**
 * Computes the partial derivative with respect to u (i.e., x) of the Berstein-Bezier function represented by the matrix
 * c. Works for arbitrary degree's (i.e., for any matrix c of arbitrary size).
 * @param c The matrix representation of the Berstein-Bezier function
 * @return The matrix representation of the row derivative of the given Berstein-Bezier function.
 * The number of row will decrease by one.
 */
Matrix row_derivative(const Matrix& c)
{
    Matrix result(c.rows() - 1, c.cols());
    for(int i = 0; i < c.rows() - 1; i++)
        for(int j = 0; j < c.cols(); j++)
            result(i, j) = (c.rows() - 1) * (c(i + 1, j) - c(i, j));
    return result;
}

/**
 * Computes the partial derivative with respect to v (i.e., y) of the Berstein-Bezier function represented by the matrix
 * c. Works for arbitrary degree's (i.e., for any matrix c of arbitrary size).
 * @param c The matrix representation of the Berstein-Bezier function
 * @return The matrix representation of the column derivative of the given Berstein-Bezier function.
 * The number of columns will decrease by one.
 */
Matrix col_derivative(const Matrix& c)
{
    Matrix result(c.rows(), c.cols() - 1);
    for(int i = 0; i < c.rows(); i++)
        for(int j = 0; j < c.cols() - 1; j++)
            result(i, j) = (c.cols() - 1) * (c(i, j + 1) - c(i, j));
    return result;
}

/**
 * Evaluates the Berstein-Bezier function represented by the Matrix c at all Gauss Points.
 * Gauss Points consist of each of the 25 pairs of points from the set of 5 Gauss Points on the domain [0, 1].
 * Works for arbitrary degree's (i.e., for any matrix c of arbitrary size).
 *
 * @param c The matrix representation of the Berstein-Bezier function to be evaluated
 * @return A 5x5 matrix representing the evaluations of the Berstein-Bezier function at each pair of Guass Points
 * (the row represents the x Gauss Point and the column represents the y Gauss Point).
 */
Matrix evaluate(const Matrix& c)
{
    int num_gauss = gauss_points.size();
    Matrix result(num_gauss, num_gauss);
    for(int i = 0; i < num_gauss; i++)
        for(int j = 0; j < num_gauss; j++)
            result(i, j) = evaluate(c, gauss_points[i], gauss_points[j]);
    return result;
}

/**
 * Evaluates the Berstein-Bezier function represented by the Matrix c at the point (u, v).
 * Works for arbitrary degree's (i.e., for matrices of any size).
 *
 * @param c The matrix representation of the Berstein-Bezier function to be evaluated
 * @param u The x-value (must be within [0, 1]) of the point to be evaluated
 * @param v The y-value (must be within [0, 1]) of the point to be evaluated
 * @return The value of the evaluation at the chosen point
 */
double evaluate(const Matrix& c, double u, double v)
{
    int deg_left = c.rows() - 1;
    int deg_right = c.cols() - 1;
    Matrix b_left(1, deg_left + 1);
    Matrix b_right(deg_right + 1, 1);
    for (int i = 0; i < deg_left + 1; i++)
        b_left(0, i) = b_func(deg_left, i, u);
    for (int i = 0; i < deg_right + 1; i++)
        b_right(i, 0) = b_func(deg_right, i, v);
    return ((b_left * c) * b_right)(0, 0);
}

/**
 * Evaluates a real-valued function on the unit square at all Gauss Points.
 * Gauss Points consist of each of the 25 pairs of points from the set of 5 Gauss Points on the domain [0, 1].
 *
 * @param func The function f: [0, 1]^2 -> R to be evaluated
 * @return A 5x5 matrix representing the evaluations of function at each of the Guass Points
 */
MatrixXd evaluate_functional(std::function<double(double, double)> func)
{
    int num_gauss = gauss_points.size();
    MatrixXd values = MatrixXd::Zero(5, 5);
    for(int i = 0; i < num_gauss; i++){
        for(int j = 0; j < num_gauss; j++){
            values(i, j) = func(gauss_points[i], gauss_points[j]);
        }
    }
    return values;
}

/**
 * Applies Gaussian quadrature to the Berstein-Bezier function represented by the Matrix c.
 * Gauss Points used consist of each of the 25 pairs of points from the set of 5 Gauss Points on the domain [0, 1].
 *
 * @param c The matrix representation of the Berstein-Bezier function to which Gaussian quadrature will be applied
 * @return The resulting value of the Guassian quadrature of the given Berstein-Bezier function.
 */
double quadrature(const Matrix& c)
{
    int num_gauss = gauss_weights.size();
    double result;
    for(int i = 0; i < num_gauss; i++)
        for(int j = 0; j < num_gauss; j++)
            result += 0.25 * gauss_weights[i] * gauss_weights[j] * c(i, j);
    return result;
}

/**
 * b function (used in evaluations)
 * @param deg degree
 * @param idx index
 * @param t evaluation point (must be within [0, 1])
 * @return
 */
double b_func(int deg, int idx, double t)
{
    return n_choose_k(deg, idx) * pow(1 - t, deg - idx) * pow(t, idx);
}

/**
 * n choose k
 * @param n must be greater than or equal to 0
 * @param k must be within {0, 1, ..., n}
 * @return returns n choose
 */
double n_choose_k(int n, int k)
{
    if(k == 0){
        return 1;
    }
    return n * n_choose_k(n - 1, k - 1) / k;
}

/**
 * Raises the degree of the matrix representation of a Berstein-Bezier function to be bi-cubic
 * @param m The matrix representation of the given Berstein-Bezier function
 * @return The bi-cubic matrix representation of the given Berstein-Bezier function
 */
Matrix degRaise(Matrix& m)
{
    int m_DegU = m.rows() - 1;
    int m_DegV = m.cols() - 1;

    bool t_IsRaiseU = m_DegU < 3 ? true : false;
    bool t_IsRaiseV = m_DegV < 3 ? true : false;

    // This patch doesn't need deg raise
    if(!t_IsRaiseU && !t_IsRaiseV)
    {
        return m;
    }

    int t_UR = t_IsRaiseU ? m_DegU+1 : m_DegU;
    int t_VR = t_IsRaiseV ? m_DegV+1 : m_DegV;

    Matrix t_BBUR(t_UR+1, m_DegV+1); // u-direction deg raised BB-coef

    // deg raise u direction
    if(t_IsRaiseU)
    {
        for(int i=0; i<=t_UR; i++)
        {
            int k = t_UR - i;
            for(int j=0; j<=m_DegV; j++)
            {
                int a = (i-1 < 0) ? 0 : i-1;
                int b = (i > m_DegU) ? i-1 : i;
                t_BBUR(i, j) = (i*m(a, j) + k*m(b, j)) / t_UR;
            }
        }
        if(!t_IsRaiseV)
        {
            return t_BBUR;
        }
    }

    // deg raise v direction
    if(t_IsRaiseV && !t_IsRaiseU)
    {
        Matrix t_BBVR(m_DegU+1, t_VR+1); // v-direction deg raised BB-coef
        for(int j=0; j<=t_VR; j++)
        {
            int k = t_VR - j;
            for(int i=0; i<=m_DegU; i++)
            {
                int a = (j-1 < 0) ? 0 : j-1;
                int b = (j > m_DegU) ? j-1 : j;
                t_BBVR(i, j) = (j*m(i, a) + k*m(i, b)) / t_VR;
            }
        }
        return t_BBVR;
    }
    else  // deg raise u & v direction
    {
        Matrix t_BBR(t_UR+1, t_VR+1); // both-direction deg raised BB-coef
        for(int j=0; j<=t_VR; j++)
        {
            int k = t_VR - j;
            for(int i=0; i<=t_UR; i++)
            {
                int a = (j-1 < 0) ? 0 : j-1;
                int b = (j > t_UR) ? j-1 : j;
                t_BBR(i, j) = (j*t_BBUR(i, a) + k*t_BBUR(i, b)) / t_VR;
            }
        }
        return t_BBR;
    }
}



// Others

template <typename T> std::vector<T> duplicate_vector(int a_Times, const std::vector<T>& a_Vector)
{
    std::vector<T> t_DuplicatedVector;
    for(int i=0; i<a_Times; i++)
    {
        t_DuplicatedVector.insert(t_DuplicatedVector.end(), a_Vector.begin(), a_Vector.end());
    }
    return t_DuplicatedVector;
}
template std::vector<int> duplicate_vector(int a_Times, const std::vector<int>& a_Vector);
template std::vector<double> duplicate_vector(int a_Times, const std::vector<double>& a_Vector);


} // end of namespace Helper
