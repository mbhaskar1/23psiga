/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PatchConstructor.hpp"
#include "../Helper/Helper.hpp"

typedef Matrix Mat48x12d;
typedef Matrix Mat80x20d;
typedef Matrix Mat384x24d;
typedef Matrix Mat448x28d;
typedef Matrix Mat512x32d;

class NGonPatchConstructor : public PatchConstructor
{
public:
    NGonPatchConstructor(const MeshType& a_Mesh)
        : m_Mesh(a_Mesh),
          m_MaskSct3(getMaskSct3()),
          m_MaskSct5(getMaskSct5()),
          m_MaskSct6(getMaskSct6()),
          m_MaskSct7(getMaskSct7()),
          m_MaskSct8(getMaskSct8()) {};

    bool isSamePatchType(const FaceHandle& a_FaceHandle) override;
    std::vector<VertexHandle> initNeighborVerts(const FaceHandle& a_FaceHandle) override;
    std::vector<Patch> getPatch(const FaceHandle& a_FaceHandle) override;
    std::vector<Patch> getPatch(const std::vector<VertexHandle>& a_NBVertexHandles) override;

    Matrix getMask() override;
    std::string getName() override { return "ngonSct" + std::to_string(m_FaceValence); };
    int getDegU() override { return 3; };
    int getDegV() override { return 3; };

private:
    const MeshType& m_Mesh;
    const Mat48x12d m_MaskSct3;
    const Mat80x20d m_MaskSct5;
    const Mat384x24d m_MaskSct6;
    const Mat448x28d m_MaskSct7;
    const Mat512x32d m_MaskSct8;
    int m_FaceValence;

    Mat48x12d getMaskSct3();
    Mat80x20d getMaskSct5();
    Mat384x24d getMaskSct6();
    Mat448x28d getMaskSct7();
    Mat512x32d getMaskSct8();
};
