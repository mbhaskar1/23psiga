/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

//  C++ std library includes
#include <iostream>
#include <string>

//  Open Mesh includes
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

//  Our own headers
#include "Pool/Pool.hpp"
#include "PatchConsumer/PatchConsumer.hpp"
#include "PatchConsumer/BVWriter.hpp"
#include "PatchConsumer/IGSWriter.hpp"
#include "ProcessMesh.hpp"

typedef OpenMesh::PolyMesh_ArrayKernelT<> MeshType;

int main(int argc, char **argv)
{
    // Check correct number of input arguments
    if(argc < 2)
    {
        std::cout << "\nArgument error!\n";
        std::cout << "Correct usage:`./PolyhedralSplines [OPTIONS] ${INPUT_FILENAME}`\n";
        std::cout << "Options: \n";
        std::cout << "-d --DEGREE_RAISE : raise deg 2 patches to deg 3. \n";
        return 1;
    }


    bool t_IsDegRaise = false;

    for (int i = 1; i < argc-1; i++)
    {
        if(std::string(argv[i]) == "--DEGREE_RAISE" || std::string(argv[i]) == "-d")
        {
            t_IsDegRaise = true;
        }
        else
        {
            std::cout << "Invalid arguments, please try again.\n";
        }
    }


    // Load mesh from .obj file
    MeshType t_Mesh;
    const std::string t_InputFile = argv[argc-1];
    OpenMesh::IO::read_mesh(t_Mesh, t_InputFile);

    // Init .bv file writer
    const std::string t_FileName = "output.bv";
    PatchConsumer* t_Writer = new BVWriter(t_FileName);

    // Users can comment out BVWriter and use IGSWriter
    // or implement MyWriter to have customized output for PDE solver, momentum computation, etc.
    /*const std::string t_FileName = "output.igs";
    PatchConsumer* t_Writer = new IGSWriter(t_FileName);*/

    // Convert mesh into Patches (contain BB-coefficients) and write patches into .bv file
    MeshInfo meshInfo = process_mesh(t_Mesh, t_IsDegRaise);

    t_Writer->start();
    for(const auto& patch : meshInfo.patches)
    {
        t_Writer->consume(patch);
    }
    t_Writer->stop();

    auto patches = meshInfo.patches;
    auto vertex_supports = meshInfo.vertex_supports;
    auto patch_supports = meshInfo.patch_supports;
    auto evaluated_bases = meshInfo.evaluated_bases;
    auto evaluated_row_diff = meshInfo.evaluated_row_diff;
    auto evaluated_col_diff = meshInfo.evaluated_col_diff;

    std::map<int, double> u_map = IGA::HeatSimulation(t_Mesh, patches,
                                                         vertex_supports, patch_supports,
                                                         evaluated_bases, evaluated_row_diff,
                                                         evaluated_col_diff, 0.5, 50);

    std::vector<Matrix> u_bb = process_mesh_values(t_Mesh, u_map, t_IsDegRaise);

    int N = 10;
    Matrix surf = Matrix::Zero((N + 1) * (N + 1) * patches.size(), 4);

    for (int patch_idx = 0; patch_idx < patches.size(); patch_idx++) {
        Matrix bb_x = patches[patch_idx].m_BBcoefs[0];
        Matrix bb_y = patches[patch_idx].m_BBcoefs[1];
        Matrix bb_z = patches[patch_idx].m_BBcoefs[2];
        Matrix bb_v = u_bb[patch_idx];

        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= N; j++) {
                int idx = patch_idx * (N + 1) * (N + 1) + i * (N + 1) + j;
                surf(idx, 0) = Helper::evaluate(bb_x, 1.0 * i / N, 1.0 * j / N);
                surf(idx, 1) = Helper::evaluate(bb_y, 1.0 * i / N, 1.0 * j / N);
                surf(idx, 2) = Helper::evaluate(bb_z, 1.0 * i / N, 1.0 * j / N);
                surf(idx, 3) = Helper::evaluate(bb_v, 1.0 * i / N, 1.0 * j / N);
            }
        }
    }

    std::cout << "Computed surface" << std::endl;
    std::cout << "Total Num Points: " << surf.rows() << std::endl;

    std::ofstream outFile;
    outFile.open("./surface_t_50.sf");
    for (int i = 0; i < patches.size() * (N + 1) * (N + 1); i++) {
        outFile << surf(i, 0) << " " << surf(i, 1) << " " << surf(i, 2) << " " << surf(i, 3) << std::endl;
    }
    outFile.close();

    delete t_Writer;

    return 0;
}
