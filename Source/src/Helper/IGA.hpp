/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
#include "../Patch/Patch.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef MatrixXd Matrix;
typedef VectorXd Vector;

typedef OpenMesh::PolyMesh_ArrayKernelT<> MeshType;
typedef MeshType::VertexHandle VertexHandle;
typedef MeshType::EdgeHandle EdgeHandle;
typedef MeshType::FaceHandle FaceHandle;
typedef MeshType::HalfedgeHandle HalfedgeHandle;
typedef MeshType::Point Point;



namespace IGA
{
    struct G_Info {
        std::vector<Matrix> G11;
        std::vector<Matrix> G12;
        std::vector<Matrix> G21;
        std::vector<Matrix> G22;
        std::vector<Matrix> det;
    };

    G_Info compute_G_alphas(std::vector<Patch> patches);

    Matrix compute_K(G_Info &evaluated_G_alphas,
                       std::map<int, std::set<int>> vertex_supports,
                       std::map<int, std::set<int>> patch_supports,
                       std::map<std::pair<int, int>, Matrix> evaluated_row_diff,
                       std::map<std::pair<int, int>, Matrix> evaluated_col_diff);

    Matrix compute_M(G_Info &evaluated_G_alphas,
                       std::map<int, std::set<int>> vertex_supports,
                       std::map<int, std::set<int>> patch_supports,
                       std::map<std::pair<int, int>, Matrix> evaluated_bases);

    Vector compute_rhs(G_Info &evaluated_G_alphas,
                         std::map<int, std::set<int>> vertex_supports,
                         std::map<int, std::set<int>> patch_supports,
                         std::map<std::pair<int, int>, Matrix> evaluated_bases,
                         std::map<int, std::function<double(double, double)>> rhs_f);

    std::vector<std::set<int>> getSupport(const std::string& patch_type);
    std::map<std::pair<int, int>, Matrix> getBases(const Matrix& mask, const std::vector<std::set<int>> &vert_supports,
                                                   int degU, int degV);
    std::map<std::pair<int, int>, Matrix> getBases(const Matrix& mask, const std::vector<std::set<int>> &vert_supports);

    std::map<int, double> HeatSimulation(const MeshType &mesh,
                                         std::vector<Patch> patches,
                                         std::map<int, std::set<int>> vertex_supports,
                                         std::map<int, std::set<int>> patch_supports,
                                         std::map<std::pair<int, int>, Matrix> evaluated_bases,
                                         std::map<std::pair<int, int>, Matrix> evaluated_row_diff,
                                         std::map<std::pair<int, int>, Matrix> evaluated_col_diff,
                                         double dest_time,
                                         int n_steps);
}
