/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "Patch.hpp"

typedef OpenMesh::PolyMesh_ArrayKernelT<> MeshType;
typedef MeshType::VertexHandle VertexHandle;
typedef MeshType::FaceHandle FaceHandle;

class PatchConstructor
{
public:
    virtual bool isSamePatchType(const VertexHandle&) { return false; };
    virtual bool isSamePatchType(const FaceHandle&) { return false; };
    virtual std::vector<VertexHandle> initNeighborVerts(const VertexHandle&) { return {}; };
    virtual std::vector<VertexHandle> initNeighborVerts(const FaceHandle&) { return {}; };
    virtual std::vector<Patch> getPatch(const VertexHandle&) { return {}; };
    virtual std::vector<Patch> getPatch(const FaceHandle&) { return {}; };
    virtual std::vector<Patch> getPatch(const std::vector<VertexHandle>&) { return {}; }
    virtual Matrix getMask() { return {}; };
    virtual std::string getName() { return ""; };
    virtual int getDegU() { return 3; };
    virtual int getDegV() { return 3; };

    std::vector<Matrix> getValues(Vector& coefficients) {
        Matrix mask = getMask();
        int degU = getDegU();
        int degV = getDegV();
        int num_coeff_per_patch = (degU + 1) * (degV + 1);

        Vector t_BBvalues = mask * coefficients;

        std::vector<Matrix> bb_values;
        for(int patch_idx = 0; patch_idx < t_BBvalues.size() / num_coeff_per_patch; patch_idx++){
            bb_values.emplace_back(degU + 1, degV + 1);
            for(int i = 0; i < degU + 1; i++){
                for(int j = 0; j < degV + 1; j++){
                    bb_values[patch_idx](i, j) =
                            t_BBvalues(patch_idx * num_coeff_per_patch + i * (degV + 1) + j, 0);
                }
            }
        }

        return bb_values;
    }
};
