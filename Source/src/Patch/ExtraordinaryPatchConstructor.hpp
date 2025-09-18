/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PatchConstructor.hpp"
#include "../Helper/Helper.hpp"

typedef Matrix Mat48x7d;
typedef Matrix Mat80x11d;
typedef Matrix Mat384x13d;
typedef Matrix Mat448x15d;
typedef Matrix Mat512x17d;

class ExtraordinaryPatchConstructor : public PatchConstructor
{
public:
    ExtraordinaryPatchConstructor(const MeshType& a_Mesh)
        : m_Mesh(a_Mesh),
          m_MaskSct3(getMaskSct3()),
          m_MaskSct5(getMaskSct5()),
          m_MaskSct6(getMaskSct6()),
          m_MaskSct7(getMaskSct7()),
          m_MaskSct8(getMaskSct8()) {};

    bool isSamePatchType(const VertexHandle& a_VertexHandle) override;
    std::vector<VertexHandle> initNeighborVerts(const VertexHandle& a_VertexHandle) override;
    std::vector<Patch> getPatch(const VertexHandle& a_VertexHandle) override;
    std::vector<Patch> getPatch(const std::vector<VertexHandle>& a_NBVertexHandles) override;

    Matrix getMask() override;
    std::string getName() override { return "eopSct" + std::to_string(m_Valence); };
    int getDegU() override { return 3; };
    int getDegV() override { return 3; };

private:
    const MeshType& m_Mesh;
    const Mat48x7d m_MaskSct3;
    const Mat80x11d m_MaskSct5;
    const Mat384x13d m_MaskSct6;
    const Mat448x15d m_MaskSct7;
    const Mat512x17d m_MaskSct8;
    int m_Valence = 0;

    Mat48x7d getMaskSct3();
    Mat80x11d getMaskSct5();
    Mat384x13d getMaskSct6();
    Mat448x15d getMaskSct7();
    Mat512x17d getMaskSct8();
};
