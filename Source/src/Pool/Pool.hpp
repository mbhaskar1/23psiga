/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "../Patch/Patch.hpp"
#include "../Patch/T0PatchConstructor.hpp"
#include "../Patch/T1PatchConstructor.hpp"
#include "../Patch/T2PatchConstructor.hpp"
#include "../Patch/ExtraordinaryPatchConstructor.hpp"
// #include "../Patch/TwoTrianglesTwoQuadsPatchConstructor.hpp"
#include "../Patch/NGonPatchConstructor.hpp"
#include "../Patch/PolarPatchConstructor.hpp"
#include "../Patch/RegularPatchConstructor.hpp"

/*
 * Store the pointers of patch constructor
 */
class PatchConstructorPool
{
public:
    PatchConstructorPool(const MeshType& a_Mesh)
    {
        // Functions to generate patches
        m_PatchConstructorPool.push_back(new RegularPatchConstructor{a_Mesh});
        m_PatchConstructorPool.push_back(new ExtraordinaryPatchConstructor{a_Mesh});
        // m_PatchConstructorPool.push_back(new TwoTrianglesTwoQuadsPatchConstructor{a_Mesh}); // deprecated
        m_PatchConstructorPool.push_back(new PolarPatchConstructor{a_Mesh});
        m_PatchConstructorPool.push_back(new T0PatchConstructor{a_Mesh});
        m_PatchConstructorPool.push_back(new T1PatchConstructor{a_Mesh});
        m_PatchConstructorPool.push_back(new T2PatchConstructor{a_Mesh});
        m_PatchConstructorPool.push_back(new NGonPatchConstructor{a_Mesh});
    }
    template <typename T> PatchConstructor* getPatchConstructor(const T& a_T) const;

private:
    std::vector<PatchConstructor*> m_PatchConstructorPool;
};

template <typename T> PatchConstructor* PatchConstructorPool::getPatchConstructor(const T& a_T) const
{
    for(auto t_PatchConstructor : m_PatchConstructorPool)
    {
        if(t_PatchConstructor->isSamePatchType(a_T))
        {
            return t_PatchConstructor;
        }
    }
    return nullptr;
}
template PatchConstructor* PatchConstructorPool::getPatchConstructor(const VertexHandle& a_VertexHandle) const;
template PatchConstructor* PatchConstructorPool::getPatchConstructor(const FaceHandle& a_FaceHandle) const;
