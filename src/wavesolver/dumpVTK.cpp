
#include <vector>
#include <map>
#include <sstream>
#include <Eigen/Dense>

#include "bubbles/FileInput.hpp"
#include "bubbles/Mesh.hpp"

#include <vtk/vtkSmartPointer.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkUnstructuredGrid.h>
#include <vtk/vtkDoubleArray.h>
#include <vtk/vtkCellType.h>
#include <vtk/vtkCellData.h>
#include <vtk/vtkXMLUnstructuredGridWriter.h>
#include <vtk/vtkTriangle.h>

using namespace std;

typedef Eigen::Matrix<double, 1, 1> VelocityValue;
typedef std::map<int, std::vector<VelocityValue>> SurfaceVelocityData;

void writeVTKFile(double t, const FileNames& files)
{
    Mesh m;
    m.loadGmsh(files.meshFile);
    std::vector<BubbleInputInfo> b = parseFreqFile(files.freqFile);

    SurfaceVelocityData v = loadSurfaceDatFile(b,
                                               files.datFile,
                                               m);


    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (auto& p : m.m_vertices)
    {
        points->InsertNextPoint(p(0), p(1), p(2));
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);

    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

    vtkIdType ids[3];
    for (int i = 0; i < m.m_triangles.size(); ++i)
    {
        auto& t = m.m_triangles[i];
        ids[0] = t(0);
        ids[1] = t(1);
        ids[2] = t(2);

        unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, ids);
    }

    vector<vtkSmartPointer<vtkDoubleArray>> cellArrays;
    for (auto& vb : v)
    {
        cellArrays.push_back(vtkSmartPointer<vtkDoubleArray>::New());
        vtkSmartPointer<vtkDoubleArray> vel = cellArrays.back();

        ostringstream os;
        os << "vel_" << vb.first;

        vel->SetName(os.str().c_str());

        //int ind = 0;

        for (int i = 0; i < m.m_triangles.size(); ++i)
        {
            if (m.m_triType[i] == Mesh::FLUID_AIR)
            {
                vel->InsertNextValue(vb.second.at(m.m_fullToSurf[i])(0));
                //vel->InsertNextValue(vb.second.at(ind)(0));
                //ind++;
            }
            else
            {
                vel->InsertNextValue(0);
            }
        }

        unstructuredGrid->GetCellData()->AddArray(vel);
    }

    ostringstream os;
    os << "data_" << t << ".vtu";

    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(os.str().c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->Write();
}

int main(int argc, char** argv)
{
    std::map<double, FileNames> files = parseFileNames(argv[1]);

    for (auto& k : files)
    {
        std::cout << k.second.meshFile << std::endl;
        writeVTKFile(k.first, k.second);
    }

    return 0;
}

