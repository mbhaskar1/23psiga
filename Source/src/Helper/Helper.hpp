/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>

using Eigen::MatrixXd;

typedef MatrixXd Matrix;

typedef OpenMesh::PolyMesh_ArrayKernelT<> MeshType;
typedef MeshType::VertexHandle VertexHandle;
typedef MeshType::EdgeHandle EdgeHandle;
typedef MeshType::FaceHandle FaceHandle;
typedef MeshType::HalfedgeHandle HalfedgeHandle;
typedef MeshType::Point Point;

typedef std::vector<double> Vec3d;



namespace Helper
{

// Vert functions
int get_vert_valence(const MeshType& a_Mesh, const VertexHandle& a_VertexHandle);
bool is_vert_3_valence(const MeshType& a_Mesh, const VertexHandle& a_VertexHandle);
bool is_vert_4_valence(const MeshType& a_Mesh, const VertexHandle& a_VertexHandle);
bool is_vert_5_valence(const MeshType& a_Mesh, const VertexHandle& a_VertexHandle);
bool are_verts_of_face_all_4_valence(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle);
void set_vert_vector_to_default(const int a_Size, std::vector<VertexHandle>& a_VertexHandles);
std::vector<FaceHandle> get_faces_around_vert_counterclock(const MeshType& a_Mesh, const VertexHandle& a_VertHandle);
std::vector<FaceHandle> get_two_layers_faces_around_vert(const MeshType& a_Mesh, const VertexHandle& a_VertHandle);
std::vector<VertexHandle> get_two_layers_verts_around_vert(const MeshType& a_Mesh, const VertexHandle& a_VertHandle);

// Face functions
bool is_triangle(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle);
bool is_quad(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle);
bool is_pentagon(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle);
bool is_hexagon(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle);
std::vector<FaceHandle> init_neighbor_faces(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle);
bool has_7_neighbor_faces(const std::vector<FaceHandle>& a_NBFaceHandles);
bool has_9_neighbor_faces(const std::vector<FaceHandle>& a_NBFaceHandles);
bool are_faces_all_quads(const MeshType& a_Mesh, const std::vector<FaceHandle>& a_FaceHandles);
int get_num_of_verts_for_face(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle);
int get_num_of_neighbor_faces(const std::vector<FaceHandle>& a_NBFaceHandles);
std::vector<VertexHandle> get_verts_of_face(const MeshType& a_Mesh, const FaceHandle& a_FaceHandle);
int num_of_quads(const MeshType& a_Mesh, std::vector<FaceHandle> a_FaceHandles);
int num_of_triangles(const MeshType& a_Mesh, std::vector<FaceHandle> a_FaceHandles);

// Type conversion
Vec3d verthandles_to_point_vec(const MeshType& a_Mesh, const VertexHandle& a_VertHandle);
Matrix verthandles_to_points_mat(const MeshType& a_Mesh, const std::vector<VertexHandle>& a_VertHandle);

// Gauss Point Evaluation
// Gauss points scaled and shifted onto [0, 1]
const std::vector<double> gauss_points(
        {0.5 + 0.5 * -0.9061798459386640,
         0.5 + 0.5 * -0.5384693101056831,
         0.5 + 0.5 * 0.0,
         0.5 + 0.5 * 0.5384693101056831,
         0.5 + 0.5 * 0.9061798459386640});
const std::vector<double> gauss_weights(
        {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891});

Matrix row_derivative(const Matrix &c);
Matrix col_derivative(const Matrix &c);
Matrix evaluate(const Matrix &c);
double evaluate(const Matrix &c, double u, double v);
Matrix evaluate_functional(std::function<double(double, double)> func);
double quadrature(const Matrix &c);
double b_func(int deg, int idx, double t);
double n_choose_k(int n, int k);
Matrix degRaise(Matrix& m);

// Others
template <typename T> std::vector<T> duplicate_vector(int a_Times, const std::vector<T>& a_Vector);

} // end of Helper namespace
