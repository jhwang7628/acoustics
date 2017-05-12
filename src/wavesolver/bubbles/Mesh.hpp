#ifndef BUBBLE_MESH_HPP
#define BUBBLE_MESH_HPP

#include <string>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

// Class to hold a Gmsh triangle mesh
class Mesh
{
public:
    enum
    {
        FLUID_AIR = 1,
        SOLID = 2
    };

	std::vector<Eigen::Vector3d> m_vertices;
	std::vector<Eigen::Vector3i> m_triangles;
    std::vector<Eigen::Vector3d> m_allTriCenters;
    std::vector<int> m_triType;

    std::vector<int> m_surfTris; // list of surface triangles (where the velocity solution data is)
    std::vector<Eigen::Vector3d> m_surfTriCenters;
    std::map<int, int> m_fullToSurf;

    void
	loadGmsh(const std::string &fileName)
	{
	    m_vertices.clear();
	    m_triangles.clear();
	    m_triType.clear();
	    m_surfTris.clear();
	    m_surfTriCenters.clear();

        std::ifstream in(fileName.c_str());

		std::string line;

		// Skip first four lines
		for (int i = 0; i < 4; ++i)
		{
			std::getline(in, line);
		}

		// Next line is # of vertices
		int numVerts;
		in >> numVerts;

		m_vertices.resize(numVerts);

		// Read the vertices
		for (int i = 0; i < numVerts; ++i)
		{
		    int index;
		    double x, y, z;

		    in >> index >> x >> y >> z;

            m_vertices[i] << x, y, z;
		}

		// Skip two lines
		in >> line >> line;

		// Read # of triangles
		int numFaces;
		in >> numFaces;

		m_triangles.resize(numFaces);
		m_triType.resize(numFaces);

		// Read triangles
		for (int i = 0; i < numFaces; ++i)
        {
            int ignore, type, index, v1, v2, v3;
            in >> ignore >> ignore >> ignore;
            in >> type >> index >> v1 >> v2 >> v3;

            m_triangles[i] << v1 - 1, v2 - 1, v3 - 1;
            m_triType[i] = type;

            // TODO: confirm whether solution data is full mesh
            // or just fluid surface data
            if (type == FLUID_AIR || type == SOLID)
            {
                m_surfTris.push_back(i);

                m_surfTriCenters.push_back( 1./3. * (m_vertices[v1] + m_vertices[v2] + m_vertices[v3]) );

                m_fullToSurf[i] = m_surfTris.size() - 1;
            }
            else
            {
                m_fullToSurf[i] = -1;
            }

            m_allTriCenters.push_back( 1./3. * (m_vertices[v1] + m_vertices[v2] + m_vertices[v3]) );
        }
	}
};

#endif // BUBBLE_MESH_HPP

