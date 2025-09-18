/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PatchConstructor.hpp"
#include "../Helper/Helper.hpp"

typedef Matrix Mat9x9d;
typedef Matrix Mat9x3d;

class TwoTrianglesTwoQuadsPatchConstructor : public PatchConstructor
{
public:
    TwoTrianglesTwoQuadsPatchConstructor(const MeshType& a_Mesh) : m_Mesh(a_Mesh) {};

    bool isSamePatchType(const VertexHandle& a_VertexHandle) override;
    std::vector<VertexHandle> initNeighborVerts(const VertexHandle& a_VertexHandle) override;
    std::vector<Patch> getPatch(const VertexHandle& a_VertexHandle) override;
    std::vector<Patch> getPatch(const std::vector<VertexHandle>& a_NBVertexHandles) override;

    Matrix getMask() override;
    std::string getName() override { return "TwoTrianglesTwoQuads"; };
    int getDegU() override { return 2; };
    int getDegV() override { return 2; };

private:
    const MeshType& m_Mesh;
    const Mat9x9d m_Mask = Matrix({
        {0.25, 0.25, 0, 0.25, 0.25, 0, 0, 0, 0},
        {0, 0.5, 0, 0, 0.5, 0, 0, 0, 0},
        {0, 0.25, 0.25, 0, 0.25, 0.25, 0, 0, 0},
        {0, 0, 0, 0.5, 0.5, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0.5, 0.5, 0, 0, 0},
        {0, 0, 0, 0.25, 0.25, 0, 0.25, 0.25, 0},
        {0, 0, 0, 0, 0.5, 0, 0, 0.5, 0},
        {0, 0, 0, 0, 0.25, 0.25, 0, 0.25, 0.25}
    });

};
