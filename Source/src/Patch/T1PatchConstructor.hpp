/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PatchConstructor.hpp"
#include "../Helper/Helper.hpp"

typedef Matrix Mat128x18d;

class T1PatchConstructor : public PatchConstructor
{
public:
    T1PatchConstructor(const MeshType& a_Mesh)
        : m_Mesh(a_Mesh), m_Mask(getMask()) {};

    bool isSamePatchType(const FaceHandle& a_FaceHandle) override;
    std::vector<VertexHandle> initNeighborVerts(const FaceHandle& a_FaceHandle) override;
    std::vector<Patch> getPatch(const FaceHandle& a_FaceHandle) override;
    std::vector<Patch> getPatch(const std::vector<VertexHandle>& a_NBVertexHandles) override;

    Mat128x18d getMask() override;
    std::string getName() override { return "T1"; };
    int getDegU() override { return 3; };
    int getDegV() override { return 3; };

private:
    const MeshType& m_Mesh;
    const Mat128x18d m_Mask;

    bool isPentagonTjunction(const FaceHandle& a_FaceHandle);

};
