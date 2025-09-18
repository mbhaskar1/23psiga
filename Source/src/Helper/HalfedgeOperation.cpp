/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */ 

#include "HalfedgeOperation.hpp"
#include "Helper.hpp"

namespace HalfedgeOperation
{

/*
 * Halfedge operation to get one vert
 * NOTE: This function changes the input halfedge
 */
VertexHandle get_vert(const MeshType& a_Mesh, HalfedgeHandle& a_HalfedgeHandle, const std::vector<int> a_Commands)
{
    VertexHandle t_VertexHandle;
    for(auto CIt=a_Commands.begin(); CIt!=a_Commands.end(); ++CIt)
    {
        switch(*CIt){
            case 1:
                a_HalfedgeHandle = a_Mesh.next_halfedge_handle(a_HalfedgeHandle);
                break;
            case 2:
                a_HalfedgeHandle = a_Mesh.prev_halfedge_handle(a_HalfedgeHandle);
                break;
            case 3:
                a_HalfedgeHandle = a_Mesh.opposite_halfedge_handle(a_HalfedgeHandle);
                break;
            case 4:
                t_VertexHandle = a_Mesh.from_vertex_handle(a_HalfedgeHandle);
                break;
        }
    }
    return t_VertexHandle;
}

/*
 * Halfedge operation to get one vert
 * NOTE: This function will not change the input halfedge
 */
VertexHandle get_vert_fixed_halfedge(const MeshType& a_Mesh, HalfedgeHandle a_HalfedgeHandle, const std::vector<int> a_Commands)
{
    return get_vert(a_Mesh, a_HalfedgeHandle, a_Commands);
}

/*
 * Operation for one corner. Operations: 1=next 2=previous 3=opposite 4=get
 * return all vertexhandles being found during the operations
 * NOTE: This function will change the input halfedge.
 */
std::vector<VertexHandle> get_verts(const MeshType& a_Mesh, HalfedgeHandle& a_HalfedgeHandle, const std::vector<int> a_Commands)
{
    std::vector<VertexHandle> t_VertexHandles;
    for(auto CIt=a_Commands.begin(); CIt!=a_Commands.end(); ++CIt)
    {
        switch(*CIt){
            case 1:
                a_HalfedgeHandle = a_Mesh.next_halfedge_handle(a_HalfedgeHandle);
                break;
            case 2:
                a_HalfedgeHandle = a_Mesh.prev_halfedge_handle(a_HalfedgeHandle);
                break;
            case 3:
                a_HalfedgeHandle = a_Mesh.opposite_halfedge_handle(a_HalfedgeHandle);
                break;
            case 4:
                t_VertexHandles.push_back(a_Mesh.from_vertex_handle(a_HalfedgeHandle));
        }
    }
    return t_VertexHandles;
}

/*
 * Halfedge operation to get verts
 * return all vertexhandles being found during the operations
 * NOTE: This function will not changle the input halfedge
 */
std::vector<VertexHandle> get_verts_fixed_halfedge(const MeshType& a_Mesh, HalfedgeHandle a_HalfedgeHandle, const std::vector<int> a_Commands)
{
    return get_verts(a_Mesh, a_HalfedgeHandle, a_Commands);
}

/*
 * Halfedge operation to get one face
 * NOTE: This function changes the input halfedge
 */
FaceHandle get_face(const MeshType& a_Mesh, HalfedgeHandle& a_HalfedgeHandle, const std::vector<int> a_Commands)
{
    FaceHandle t_FaceHandle;
    for(auto CIt=a_Commands.begin(); CIt!=a_Commands.end(); ++CIt)
    {
        switch(*CIt){
            case 1:
                a_HalfedgeHandle = a_Mesh.next_halfedge_handle(a_HalfedgeHandle);
                break;
            case 2:
                a_HalfedgeHandle = a_Mesh.prev_halfedge_handle(a_HalfedgeHandle);
                break;
            case 3:
                a_HalfedgeHandle = a_Mesh.opposite_halfedge_handle(a_HalfedgeHandle);
                break;
            case 5:
                t_FaceHandle = a_Mesh.face_handle(a_HalfedgeHandle);
                break;
        }
    }
    return t_FaceHandle;
}


/*
 * Fill the TargetVertexHandles vector by given the index ordering and commands
 */
void init_verts(const MeshType& a_Mesh, HalfedgeHandle a_HalfedgeHandle, std::vector<VertexHandle>& a_TargetVertexHandles,
                const std::vector<std::vector<int>>& a_OrderOfVertIdx, const std::vector<int>& a_CommandForOneSet)
{
    for(auto vertIdSetIt=a_OrderOfVertIdx.begin(); vertIdSetIt!=a_OrderOfVertIdx.end(); ++vertIdSetIt)
    {
        auto t_SetOfVectHandles = get_verts(a_Mesh, a_HalfedgeHandle, a_CommandForOneSet);
        for(int i=0; i<t_SetOfVectHandles.size(); i++)
        {
            a_TargetVertexHandles[vertIdSetIt->at(i)] = t_SetOfVectHandles[i];
        }
    }
}

/*
 * Fill the TargetVertexHandles vector by given the index ordering and commands
 * NOTE: This function will changle the input halfedge!
 */
void init_verts(const MeshType& a_Mesh, HalfedgeHandle& a_HalfedgeHandle, std::vector<VertexHandle>& a_TargetVertexHandles,
                const std::vector<int>& a_OrderOfVertIdx, const std::vector<int>& a_Commands)
{
    auto t_VertexHandles = get_verts(a_Mesh, a_HalfedgeHandle, a_Commands);
    for(int i=0; i<t_VertexHandles.size(); i++)
    {
        a_TargetVertexHandles[a_OrderOfVertIdx[i]] = t_VertexHandles[i];
    }
}

/*
 * Fill the TargetVertexHandles vector by given the index ordering and commands
 * NOTE: This function will changle the input halfedge!
 */
void init_verts_fixed_halfedge(const MeshType& a_Mesh, HalfedgeHandle a_HalfedgeHandle, std::vector<VertexHandle>& a_TargetVertexHandles,
                               const std::vector<int>& a_OrderOfVertIdx, const std::vector<int>& a_Commands)
{
    auto t_VertexHandles = get_verts(a_Mesh, a_HalfedgeHandle, a_Commands);
    for(int i=0; i<t_VertexHandles.size(); i++)
    {
        a_TargetVertexHandles[a_OrderOfVertIdx[i]] = t_VertexHandles[i];
    }
}

/*
 * Duplicate commands by specific times
 */
std::vector<int> duplicate_commands(int a_Times, const std::vector<int>& a_Commands)
{
    return Helper::duplicate_vector(a_Times, a_Commands);
}

} // end of HalfedgeOperation namespace
