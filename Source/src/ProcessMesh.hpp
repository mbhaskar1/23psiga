/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include "Pool/Pool.hpp"
#include "PatchConsumer/PatchConsumer.hpp"
#include "Helper/Helper.hpp"
#include "Helper/IGA.hpp"

typedef OpenMesh::PolyMesh_ArrayKernelT<> MeshType;

struct MeshInfo {
    std::vector<Patch> patches;   // Converted to bv file, i.e., vector of 4x4 matrices, which you then use to create G_alpha's
    std::map<int, std::set<int>> vertex_supports;  // List of patch indices in support for each control point
    std::map<int, std::set<int>> patch_supports;    // List of vertex indices whose support contains given patch index
    std::map<std::pair<int, int>, Matrix> evaluated_bases;  // For each vertex index and patch index (in support), evaluated Nia
    std::map<std::pair<int, int>, Matrix> evaluated_row_diff; // " ", evaluated row-derivative of Nia
    std::map<std::pair<int, int>, Matrix> evaluated_col_diff; // " ", evaluated col-derivative of Nia
};

MeshInfo process_mesh(const MeshType& a_Mesh, const bool a_IsDegRaise)
{
	// Construct the pool which will process the mesh
	PatchConstructorPool t_PatchConstructorPool(a_Mesh);

	MeshInfo meshInfo;

    // Vert iteration
	MeshType::VertexIter t_VertIt, t_VertEnd(a_Mesh.vertices_end());
	for (auto t_VertIt = a_Mesh.vertices_begin(); t_VertIt != t_VertEnd; ++t_VertIt)
	{
		auto* t_Constructor = t_PatchConstructorPool.getPatchConstructor(*t_VertIt);
		if (t_Constructor == nullptr)
		{
			continue;
		}
        std::cout << t_Constructor->getName() << std::endl;

        // TODO: THIS PREVIOUSLY CAUSED US ISSUES WHICH WERE ONLY RESOLVED WHEN WE MADE EVERYTHING BI-CUBIC
        int degU = t_Constructor->getDegU();
        int degV = t_Constructor->getDegV();

        auto t_NBVerts = t_Constructor->initNeighborVerts(*t_VertIt);
        auto vert_supports = IGA::getSupport(t_Constructor->getName());
        auto bases = IGA::getBases(t_Constructor->getMask(), vert_supports, degU, degV);
		auto t_VertPatches = t_Constructor->getPatch(t_NBVerts);

        int starting_patch_idx = meshInfo.patches.size();

        for (int patch_idx = 0; patch_idx < t_VertPatches.size(); patch_idx++) {
            meshInfo.patch_supports[starting_patch_idx + patch_idx] = std::set<int>();
        }

        for (int i = 0; i < t_NBVerts.size(); i++) {
            int vert_idx = t_NBVerts[i].idx();

            for (auto patch_idx : vert_supports[i]) {
                int shifted_patch_idx = starting_patch_idx + patch_idx;
                meshInfo.vertex_supports[vert_idx].insert(shifted_patch_idx);
                meshInfo.patch_supports[shifted_patch_idx].insert(vert_idx);

                auto basis_element = bases[{i, patch_idx}];
                if(a_IsDegRaise){
                    basis_element = Helper::degRaise(basis_element);
                }
                meshInfo.evaluated_bases[{vert_idx, shifted_patch_idx}] =
                        Helper::evaluate(basis_element);
                meshInfo.evaluated_row_diff[{vert_idx, shifted_patch_idx}] =
                        Helper::evaluate(Helper::row_derivative(basis_element));
                meshInfo.evaluated_col_diff[{vert_idx, shifted_patch_idx}] =
                        Helper::evaluate(Helper::col_derivative(basis_element));
            }
        }

        for (auto t_Patch : t_VertPatches)
		{
			if(a_IsDegRaise)
			{
				t_Patch.degRaise();
			}
            meshInfo.patches.push_back(t_Patch);
		}
	}

	// Face iteration
	MeshType::FaceIter t_FaceIt, t_FaceEnd(a_Mesh.faces_end());
	for (auto t_FaceIt = a_Mesh.faces_begin(); t_FaceIt != t_FaceEnd; ++t_FaceIt)
	{
		auto* t_Constructor = t_PatchConstructorPool.getPatchConstructor(*t_FaceIt);
		if (t_Constructor == nullptr)
		{
			continue;
		}
        std::cout << t_Constructor->getName() << std::endl;

        // TODO: THIS PREVIOUSLY CAUSED US ISSUES WHICH WERE ONLY RESOLVED WHEN WE MADE EVERYTHING BI-CUBIC
        int degU = t_Constructor->getDegU();
        int degV = t_Constructor->getDegV();

        auto t_NBVerts = t_Constructor->initNeighborVerts(*t_FaceIt);
        auto vert_supports = IGA::getSupport(t_Constructor->getName());
        auto bases = IGA::getBases(t_Constructor->getMask(), vert_supports, degU, degV);
		auto t_FacePatches = t_Constructor->getPatch(t_NBVerts);

        int starting_patch_idx = meshInfo.patches.size();

        for (int patch_idx = 0; patch_idx < t_FacePatches.size(); patch_idx++) {
            meshInfo.patch_supports[starting_patch_idx + patch_idx] = std::set<int>();
        }

        for (int i = 0; i < t_NBVerts.size(); i++) {
            int vert_idx = t_NBVerts[i].idx();

            for (auto patch_idx : vert_supports[i]) {
                int shifted_patch_idx = starting_patch_idx + patch_idx;
                meshInfo.vertex_supports[vert_idx].insert(shifted_patch_idx);
                meshInfo.patch_supports[shifted_patch_idx].insert(vert_idx);

                auto basis_element = bases[{i, patch_idx}];
                if(a_IsDegRaise){
                    basis_element = Helper::degRaise(basis_element);
                }
                meshInfo.evaluated_bases[{vert_idx, shifted_patch_idx}] =
                        Helper::evaluate(basis_element);
                meshInfo.evaluated_row_diff[{vert_idx, shifted_patch_idx}] =
                        Helper::evaluate(Helper::row_derivative(basis_element));
                meshInfo.evaluated_col_diff[{vert_idx, shifted_patch_idx}] =
                        Helper::evaluate(Helper::col_derivative(basis_element));
            }
        }

        for (auto t_Patch : t_FacePatches)
		{
			if(a_IsDegRaise)
			{
				t_Patch.degRaise();
			}
            meshInfo.patches.push_back(t_Patch);
		}
	}

	return meshInfo;
}

std::vector<Matrix> process_mesh_values(const MeshType &a_Mesh, std::map<int, double> coefficients, const bool a_IsDegRaise) {
    // Construct the pool which will process the mesh
    PatchConstructorPool t_PatchConstructorPool(a_Mesh);
    std::vector<Matrix> values;

    // Vert iteration
    MeshType::VertexIter t_VertIt, t_VertEnd(a_Mesh.vertices_end());
    // TODO: Check that this is looping in the same order as before
    for (auto t_VertIt = a_Mesh.vertices_begin(); t_VertIt != t_VertEnd; ++t_VertIt) {
        auto *t_Constructor = t_PatchConstructorPool.getPatchConstructor(*t_VertIt);
        if (t_Constructor == nullptr) {
            continue;
        }

        auto t_NBVerts = t_Constructor->initNeighborVerts(*t_VertIt);
        Vector t_NBValues(t_NBVerts.size());
        for (int i = 0; i < t_NBVerts.size(); i++) {
            t_NBValues(i) = coefficients[t_NBVerts[i].idx()];
        }
        auto t_Values = t_Constructor->getValues(t_NBValues);

        if(a_IsDegRaise){
            for (auto & t_Value : t_Values) {
                t_Value = Helper::degRaise(t_Value);
            }
        }

        values.insert(values.end(), t_Values.begin(), t_Values.end());
    }

    // Face iteration
    MeshType::FaceIter t_FaceIt, t_FaceEnd(a_Mesh.faces_end());
    for (auto t_FaceIt = a_Mesh.faces_begin(); t_FaceIt != t_FaceEnd; ++t_FaceIt) {
        auto *t_Constructor = t_PatchConstructorPool.getPatchConstructor(*t_FaceIt);
        if (t_Constructor == nullptr) {
            continue;
        }

        auto t_NBVerts = t_Constructor->initNeighborVerts(*t_FaceIt);
        Vector t_NBValues(t_NBVerts.size());
        for (int i = 0; i < t_NBVerts.size(); i++) {
            t_NBValues(i) = coefficients[t_NBVerts[i].idx()];
        }
        auto t_Values = t_Constructor->getValues(t_NBValues);

        if(a_IsDegRaise){
            for (auto & t_Value : t_Values) {
                t_Value = Helper::degRaise(t_Value);
            }
        }

        values.insert(values.end(), t_Values.begin(), t_Values.end());
    }

    return values;
}

