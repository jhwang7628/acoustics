#include "io/TetMeshReader.hpp"
#include "io/TglMeshReader.hpp"
#include "io/TglMeshWriter.hpp"
#include "geometry/TriangleMesh.hpp"
int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "**Usage: " << argv[0] << " <input_tet> <output_obj>\n";
        return 1;
    }
    TriangleMesh<double> mesh;
    FixVtxTetMesh<double> tetmesh;
    auto err = FV_TetMeshLoader_Double::load_mesh(argv[1], tetmesh);
    assert(err == SUCC_RETURN);
    tetmesh.extract_surface(&mesh);
    mesh.generate_normals();
    mesh.update_vertex_areas();
    MeshObjWriter::write(mesh, argv[2]);
    return 0;
}
