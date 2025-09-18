#include <fstream>
#include "IGA.hpp"

IGA::G_Info IGA::compute_G_alphas(std::vector<Patch> patches) {
    std::vector<Matrix> G11(patches.size());
    std::vector<Matrix> G12(patches.size());
    std::vector<Matrix> G21(patches.size());
    std::vector<Matrix> G22(patches.size());
    std::vector<Matrix> det(patches.size());

    for (int patch_idx = 0; patch_idx < patches.size(); patch_idx++) {
        Matrix cx = patches[patch_idx].m_BBcoefs[0];
        Matrix cy = patches[patch_idx].m_BBcoefs[1];
        Matrix cz = patches[patch_idx].m_BBcoefs[2];

        Matrix dx_u = Helper::evaluate(Helper::row_derivative(cx));
        Matrix dy_u = Helper::evaluate(Helper::row_derivative(cy));
        Matrix dz_u = Helper::evaluate(Helper::row_derivative(cz));

        Matrix dx_v = Helper::evaluate(Helper::col_derivative(cx));
        Matrix dy_v = Helper::evaluate(Helper::col_derivative(cy));
        Matrix dz_v = Helper::evaluate(Helper::col_derivative(cz));

        // First fundamental form (Think of it in the form of a 2x2 metric tensor) Also g21 is equivalent to g12, so it's not listed.
        // Summation of inner products of partial derivatives
        G11[patch_idx] =
                (dx_u.array() * dx_u.array()) + (dy_u.array() * dy_u.array()) + (dz_u.array() * dz_u.array());
        G12[patch_idx] =
                (dx_u.array() * dx_v.array()) + (dy_u.array() * dy_v.array()) + (dz_u.array() * dz_v.array());
        G21[patch_idx] =
                (dx_u.array() * dx_v.array()) + (dy_u.array() * dy_v.array()) + (dz_u.array() * dz_v.array());
        G22[patch_idx] =
                (dx_v.array() * dx_v.array()) + (dy_v.array() * dy_v.array()) + (dz_v.array() * dz_v.array());

        //This is the determinant of the metric tensor for this surface/patch.
        det[patch_idx] =
                G11[patch_idx].array() * G22[patch_idx].array() - G12[patch_idx].array() * G21[patch_idx].array();
    }

    return {G11, G12, G21, G22, det};
}

Matrix IGA::compute_K(IGA::G_Info &evaluated_G_alphas,
                      std::map<int, std::set<int>> vertex_supports,
                      std::map<int, std::set<int>> patch_supports,
                      std::map<std::pair<int, int>, Matrix> evaluated_row_diff,
                      std::map<std::pair<int, int>, Matrix> evaluated_col_diff) {
    Matrix K = Matrix::Zero(vertex_supports.size(), vertex_supports.size());
    for (int patch_idx = 0; patch_idx < patch_supports.size(); patch_idx++) {
        // CAN MAKE MORE EFFICIENT BY TAKING ADVANTAGE OF K BEING SYMMETRIC
        for (int i: patch_supports[patch_idx]) {
            for (int j: patch_supports[patch_idx]) {
                Matrix &Ni_1 = evaluated_row_diff[std::make_pair(i, patch_idx)];
                Matrix &Ni_2 = evaluated_col_diff[std::make_pair(i, patch_idx)];
                Matrix &Nj_1 = evaluated_row_diff[std::make_pair(j, patch_idx)];
                Matrix &Nj_2 = evaluated_col_diff[std::make_pair(j, patch_idx)];

                Matrix &G_11 = evaluated_G_alphas.G11[patch_idx];
                Matrix &G_12 = evaluated_G_alphas.G12[patch_idx];
                Matrix &G_21 = evaluated_G_alphas.G21[patch_idx];
                Matrix &G_22 = evaluated_G_alphas.G22[patch_idx];

                Matrix &det_G = evaluated_G_alphas.det[patch_idx];

                Matrix Ni_G_1 =
                        (Ni_1.array() * G_22.array() - Ni_2.array() * G_12.array()) / det_G.array().sqrt();
                Matrix Ni_G_2 =
                        ((-1) * Ni_1.array() * G_21.array() + Ni_2.array() * G_11.array()) / det_G.array().sqrt();

                Matrix Ni_G_Nj = Ni_G_1.array() * Nj_1.array() + Ni_G_2.array() * Nj_2.array();

                K(i, j) += Helper::quadrature(Ni_G_Nj);
            }
        }
    }

    return K;
}

Matrix IGA::compute_M(G_Info &evaluated_G_alphas,
                      std::map<int, std::set<int>> vertex_supports,
                      std::map<int, std::set<int>> patch_supports,
                      std::map<std::pair<int, int>, Matrix> evaluated_bases) {
    Matrix M = Matrix::Zero(vertex_supports.size(), vertex_supports.size());
    for (int patch_idx = 0; patch_idx < patch_supports.size(); patch_idx++) {
        // CAN MAKE MORE EFFICIENT BY TAKING ADVANTAGE OF M BEING SYMMETRIC
        for (int i: patch_supports[patch_idx]) {
            for (int j: patch_supports[patch_idx]) {
                Matrix &Ni = evaluated_bases[std::make_pair(i, patch_idx)];
                Matrix &Nj = evaluated_bases[std::make_pair(j, patch_idx)];

                Matrix &det_G = evaluated_G_alphas.det[patch_idx];

                M(i, j) += Helper::quadrature(Ni.array() * Nj.array() * det_G.array().sqrt());
            }
        }
    }

    return M;
}

Vector IGA::compute_rhs(G_Info &evaluated_G_alphas,
                        std::map<int, std::set<int>> vertex_supports,
                        std::map<int, std::set<int>> patch_supports,
                        std::map<std::pair<int, int>, Matrix> evaluated_bases,
                        std::map<int, std::function<double(double, double)>> rhs_f) {
    Vector rhs = Vector::Zero(vertex_supports.size());
    for (int patch_idx = 0; patch_idx < patch_supports.size(); patch_idx++) {
        Matrix evaluated_rhs_f = Helper::evaluate_functional(rhs_f[patch_idx]);
        Matrix &det_G = evaluated_G_alphas.det[patch_idx];
        for (int i: patch_supports[patch_idx]) {
            Matrix &Ni = evaluated_bases[std::make_pair(i, patch_idx)];
            rhs(i) += Helper::quadrature(evaluated_rhs_f.array() * Ni.array() * det_G.array().sqrt());
        }
    }

    return rhs;
}

/**
 * For each vertex in a patch configuration, get its support, i.e., the patches which will be affected by changes made
 * to that vertex
 * @param patch_type the name of the patch configuration (use getName() function in PatchConstructor class)
 * @return A list of sets of integers. The index in the list denotes the vertex index within the patch configuration,
 * and the integers in the sets correspond to the patch indices within the patch configuration which are in the support
 * of the respective vertex.
 */
std::vector<std::set<int>> IGA::getSupport(const std::string &patch_type) {
    std::ifstream supp_file("./SupportInfo/" + patch_type + ".supp");

    std::vector<std::set<int>> support_vectors;
    int prev_idx = -1;
    while (supp_file) {
        int idx;
        supp_file >> idx;
        if (idx == prev_idx)
            break;
        prev_idx = idx;
        int num_patches;
        supp_file >> num_patches;
        std::set<int> support;
        for (int i = 0; i < num_patches; i++) {
            int patch_idx;
            supp_file >> patch_idx;
            support.insert(patch_idx);
        }
        support_vectors.push_back(support);
    }
    return support_vectors;
}

/**
 * Collects individual basis elements N_{i}_{\alpha} from mask file for some configuration and computed vertex supports
 * @param mask mask Matrix (use ReadCSV2Matrix on one of the csv files in the Table folder)
 * @param vert_supports computed vertex supports (use getSupport() function from above)
 * @param degU the degree in the u direction for patches in this configuration
 * @param degV the degree in the v direction for patches in this configuration
 * @return A map mapping vertex index i and patch index \alpha to the basis element N_{i}_{\alpha}
 */
std::map<std::pair<int, int>, Matrix> IGA::getBases(const Matrix &mask, const std::vector<std::set<int>> &vert_supports,
                                                    int degU, int degV) {
    const int numVertices = mask.cols();

    int deg1U = degU + 1;
    int deg1V = degV + 1;

    std::map<std::pair<int, int>, Matrix> bases;
    for (int vert_idx = 0; vert_idx < numVertices; vert_idx++) {
        for (int patch_idx: vert_supports[vert_idx]) {
            Matrix basis(deg1U, deg1V);
            for (int i = 0; i < deg1U * deg1V; i++) {
                basis(i / deg1U, i % deg1V) = mask(deg1U * deg1V * patch_idx + i, vert_idx);
            }
            bases[std::pair<int, int>(vert_idx, patch_idx)] = basis;
        }
    }

    return bases;
}

/**
 * Collects individual basis elements N_{i}_{\alpha} from mask file for some configuration and computed vertex supports
 * (assumes degree 3)
 * @param mask mask Matrix (use ReadCSV2Matrix on one of the csv files in the Table folder)
 * @param vert_supports computed vertex supports (use getSupport() function from above)
 * @return A map mapping vertex index i and patch index \alpha to the basis element N_{i}_{\alpha}
 */
std::map<std::pair<int, int>, Matrix>
IGA::getBases(const Matrix &mask, const std::vector<std::set<int>> &vert_supports) {
    return IGA::getBases(mask, vert_supports, 3, 3);
}

std::map<int, double> IGA::HeatSimulation(const MeshType &mesh,
                                          std::vector<Patch> patches,
                                          std::map<int, std::set<int>> vertex_supports,
                                          std::map<int, std::set<int>> patch_supports,
                                          std::map<std::pair<int, int>, Matrix> evaluated_bases,
                                          std::map<std::pair<int, int>, Matrix> evaluated_row_diff,
                                          std::map<std::pair<int, int>, Matrix> evaluated_col_diff,
                                          double dest_time,
                                          int n_steps) {
    auto evaluated_G_alphas = compute_G_alphas(patches);
    Matrix K = compute_K(evaluated_G_alphas, vertex_supports, patch_supports, evaluated_row_diff,
                           evaluated_col_diff);
    Matrix M = compute_M(evaluated_G_alphas, vertex_supports, patch_supports, evaluated_bases);

    double time_step = dest_time / n_steps;
    int fidx = 0;
    Vector ui = Vector::Zero(vertex_supports.size());
    ui(fidx) = 100;
    if (dest_time > 0) {
        for (int i = 1; i <= n_steps; i++) {
            Matrix A = 2 * M.array() + time_step * K.array();
            Vector b = (2 * M.array() - time_step * K.array()).matrix() * ui;
            ui = A.colPivHouseholderQr().solve(b);
            ui(fidx) = 100;
        }
    }
    std::cout << "Final ui:" << std::endl;
    std::cout << ui << std::endl;

    std::map<int, double> u_map;
    for (int i = 0; i < vertex_supports.size(); i++) {
        u_map.insert(std::make_pair(i, ui(i)));
    }
    return u_map;
}