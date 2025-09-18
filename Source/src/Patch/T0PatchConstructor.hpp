/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PatchConstructor.hpp"
#include "../Helper/Helper.hpp"

typedef Matrix Mat64x14d;

class T0PatchConstructor : public PatchConstructor
{
public:
    T0PatchConstructor(const MeshType& a_Mesh)
        : m_Mesh(a_Mesh), m_Mask(getMask()) {};

    bool isSamePatchType(const FaceHandle& a_FaceHandle) override;
    std::vector<VertexHandle> initNeighborVerts(const FaceHandle& a_FaceHandle) override;
    std::vector<Patch> getPatch(const FaceHandle& a_FaceHandle) override;
    std::vector<Patch> getPatch(const std::vector<VertexHandle>& a_NBVertexHandles) override;
    Mat64x14d getMask() override;
    std::string getName() override { return "T0"; };
    int getDegU() override { return 3; };
    int getDegV() override { return 3; };

private:
    const MeshType& m_Mesh;
    const Mat64x14d m_Mask;

};
