/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

typedef OpenMesh::PolyMesh_ArrayKernelT<> MeshType;
typedef MeshType::VertexHandle VertexHandle;
typedef MeshType::FaceHandle FaceHandle;
typedef MeshType::HalfedgeHandle HalfedgeHandle;

namespace HalfedgeOperation
{

VertexHandle get_vert(const MeshType& a_Mesh, HalfedgeHandle& a_HalfedgeHandle, const std::vector<int> a_Commands);
VertexHandle get_vert_fixed_halfedge(const MeshType& a_Mesh, HalfedgeHandle a_HalfedgeHandle, const std::vector<int> a_Commands);

FaceHandle get_face(const MeshType& a_Mesh, HalfedgeHandle& a_HalfedgeHandle, const std::vector<int> a_Commands);

std::vector<VertexHandle> get_verts(const MeshType& a_Mesh, HalfedgeHandle& a_HalfedgeHandle, const std::vector<int> a_Commands);
std::vector<VertexHandle> get_verts_fixed_halfedge(const MeshType& a_Mesh, HalfedgeHandle a_HalfedgeHandle, const std::vector<int> a_Commands);

void init_verts(const MeshType& a_Mesh, HalfedgeHandle a_HalfedgeHandle, std::vector<VertexHandle>& a_TargetVertexHandles,
                const std::vector<std::vector<int> >& a_OrderOfVertIdx, const std::vector<int>& a_CommandForOneSet);
void init_verts(const MeshType& a_Mesh, HalfedgeHandle& a_HalfedgeHandle, std::vector<VertexHandle>& a_TargetVertexHandles,
                const std::vector<int>& a_OrderOfVertIdx, const std::vector<int>& a_Commands);
void init_verts_fixed_halfedge(const MeshType& a_Mesh, HalfedgeHandle a_HalfedgeHandle, std::vector<VertexHandle>& a_TargetVertexHandles,
                               const std::vector<int>& a_OrderOfVertIdx, const std::vector<int>& a_Commands);

std::vector<int> duplicate_commands(int a_Times, const std::vector<int>& a_Commands);

}
