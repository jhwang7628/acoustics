#include <sstream> 
#include <fstream> 

#include <unistd.h>

#include "OpenfoamMesh.h" 
#include "IO/IO.h" 

#include "Grid.h" 
#include "STL_Utils/STL_Wrapper.h" 

#include "vtkConverter/vtkConverter.h"

/* 
 * FOLDER: 0
 * FILE: 1
 */
void CheckValidity(const string &name, int FOLDER_FILE) 
{
    switch (FOLDER_FILE) 
    {
        case 0: 
            if (!IO::ExistFolder(name))
                throw std::runtime_error("**ERROR** sought folder " + name + " does not exist"); 
            break;
        case 1: 
            if (!IO::ExistFile(name))
                throw std::runtime_error("**ERROR** sought file " + name + " does not exist"); 
            break;
    }
}

/* 
 * Make sure the class is correct, for example for compact list.
 */
void CheckFoamClass(std::ifstream &inFile, const std::string &className) 
{
    bool correctClass = false; 
    std::string line; 
    int intBuffer=0; 
    while(std::getline(inFile, line))
    {
        if (std::istringstream(line)>>intBuffer) break;
        if (IO::FindString(line,"class") && IO::FindString(line,className))
        {
            correctClass = true; 
            break;
        }
    }
    inFile.clear(); 
    inFile.seekg(0, ios::beg);

    if (!correctClass) 
        throw std::runtime_error("**ERROR** foam class "+className+" cannot be found in the file stream."); 
}

void CheckReadBracket(std::ifstream &inFile) 
{
    char bracket; 
    inFile.get(bracket); 
    if (bracket!='(' && bracket!=')' && bracket!='{' && bracket!='}')
        throw std::runtime_error("**ERROR** bracket check failed. read in char :" + std::string(&bracket));
}

void ReadChar(std::ifstream &inFile) 
{
    char charBuffer; 
    inFile.get(charBuffer); 
}

void CheckReadFirstInteger(std::ifstream &inFile, const int &compareTo)
{
    int intBuffer;
    bool foundInt = false; 
    std::string line; 
    while (std::getline(inFile,line))
    {
        if (std::istringstream(line) >> intBuffer)
        {
            foundInt = true; 
            break;
        }
    }

    if (!foundInt)
        throw runtime_error("**ERROR** cannot find integer.");

    if (intBuffer != compareTo)
        throw runtime_error("**ERROR** integer comparison mistake: " + std::to_string(intBuffer) + "->" + std::to_string(compareTo)); 
}


/* 
 * return whether the foam file is legal. the return file stream should be set
 * exactly at the data start position. 
 */
int ReadFoamHeader(std::ifstream &inFile, const std::string &objectName, const std::string zoneName="NO_ZONE")
{
    const bool isReadZone = (zoneName.compare("NO_ZONE")==0 ? false : true); 

    std::string line; 
    bool foundObject = false; 
    int N_readLines = -1; 

    while (std::getline(inFile, line))
    {
        // try to find object name after foamfile is found
        if (IO::findInitStr(line,"FoamFile") ) 
        {
            while (std::getline(inFile, line)) 
            {
                if (IO::FindString(line,objectName)) 
                {
                    foundObject = true; 
                    break;
                }
            }
        }

        std::istringstream iss(line); 
        if (iss>>N_readLines) 
        {
            char startBracket; 
            inFile.get(startBracket); 
            if ( startBracket != '(' && startBracket != '{' )
            {
                std::cerr << "**WARNING** Cannot parse foam header from file : " << objectName << ". abort\n"; 
                return -1; 
            }
            break;
        }
    }

    if (!isReadZone) 
    {
        // can terminate if no need to read zones 
        if (foundObject && N_readLines) 
            return N_readLines; 
        else 
            return -1; 
    }
    else
    {
        bool foundZone = false; 
        while (std::getline(inFile, line))
        {
            if (IO::findInitStr(line,zoneName))
            {
                foundZone = true; 
                break;
            }
        }

        while (std::getline(inFile,line))
        {
            std::istringstream iss(line); 
            if (iss>>N_readLines)
            {
                char startBracket; 
                inFile.get(startBracket); 
                if ( startBracket != '(' && startBracket != '{' )
                {
                    std::cerr << "**WARNING** Cannot parse foam header from file. abort\n"; 
                    return -1; 
                }
                break;
            }
        }

        if (foundObject && N_readLines && foundZone) 
            return N_readLines; 
        else 
            return -1; 
    }

}

inline int FlattenIndex(const int &Nx, const int &Ny, const int &Nz, const int &ind_x, const int &ind_y, const int &ind_z)
{
    return ind_x + ind_y*Nx + ind_z*Nx*Ny; 
}

///////////////////////////////////////////////////////////

void OpenfoamMesh::ReinitializeMesh() 
{
    std::cout << "reinitializing mesh by reading data from openfoam root path : " + OpenfoamMesh::_rootPath << std::endl;

    CheckValidity(OpenfoamMesh::_rootPath    ,0); 
    CheckValidity(OpenfoamMesh::_polyMeshPath,0); 

    ///// assume all read files are BINARY ///// 
    ReadFoamFile    ("points"       , 3, _pointPositions                        ); 
    ReadFoamFile_Int("cellZones"    , 1, _cellIndicies, OpenfoamMesh::_zoneName ); 
    ReadFoamFile_Int("owner"        , 1, _owners                                ); 
    ReadFoamFile_Int("neighbour"    , 1, _neighbours                            ); 
    ReadFoamFileFace("faces"        ,    _faces                                 ); 

    TrimAndCheckCellIndex();
}


int OpenfoamMesh::CountNonTriangleFaces() const
{
    int nonTriangleCount = 0;
    for (size_t ii=0; ii<_faces.Size(); ii++) 
    {
        if (_faces.N_polygons(ii)!=3)
        {
            std::cout << "N_polygon = " << _faces.N_polygons(ii) << std::endl;
            nonTriangleCount++; 
        }
    }
    return nonTriangleCount; 
}

void OpenfoamMesh::SplitNonTetrahedrons(const int &nonTetrahedronCellIndex, const std::vector<int> &facesThisCellOwns, std::vector<Tetrahedron> &splittedTetrahedrons)
{
    const int N_faces = facesThisCellOwns.size(); 
    int N_overlappedPoints = 0;
    for (int ii=0; ii<N_faces; ii++) 
        N_overlappedPoints += _faces.N_polygons(facesThisCellOwns[ii]);

    /* 
     * for triangular prism splitting, implementation follows the idea
     * of [1].
     *
     * ref [1]: JDPL Marie-Gabrielle 1999. Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra
     */ 
    if (N_faces==5 && N_overlappedPoints==18) // triangular prism case
    {
        Eigen::MatrixXi triangleFaces(2,3); 
        Eigen::MatrixXi quadFaces(3,4); 
        std::vector<bool> markedQuad(3,false);  // mark the two quads that contain the smallest index, i1

        int count_tri  = 0; 
        int count_quad = 0; 
        for (int ii=0; ii<N_faces; ii++) 
        {
            if (_faces.N_polygons(facesThisCellOwns[ii]) == 3) 
            {
                triangleFaces.row(count_tri) = _faces[facesThisCellOwns[ii]]; 
                count_tri++; 
            }
            else if (_faces.N_polygons(facesThisCellOwns[ii]) == 4)
            {
                quadFaces.row(count_quad) = _faces[facesThisCellOwns[ii]]; 
                count_quad++; 
            }
        }

        int i1 = triangleFaces.minCoeff(); 
        int i5=-1, i6=-1; 
        // find i5, i6, the two indicies across from i1 in the two quads
        for (int ii=0; ii<3; ii++) 
            for (int jj=0; jj<4; jj++) 
                if (quadFaces(ii,jj) == i1 && i5==-1) // first fill in i5 and mark the quad
                {
                    i5 = quadFaces(ii,(jj+2)%4); 
                    markedQuad[ii] = true; 
                } 
                else if (quadFaces(ii,jj) == i1 && i6==-1)  // then fill in i6 and mark the quad
                {
                    i6 = quadFaces(ii,(jj+2)%4); 
                    markedQuad[ii] = true; 
                } 

        int triangleContainsi5=-1;
        for (int ii=0; ii<2; ii++) 
            for (int jj=0; jj<3; jj++) 
            {
                if (i5 == triangleFaces(ii,jj)) 
                    triangleContainsi5 = ii; 
            }

        int i4 = triangleFaces.row(triangleContainsi5).sum() - i5 - i6;  // 4,5,6 on the same triangle. 1456 should be a tet

        // TODO can do more work here to figure out the min index in the last
        // quad face and use it to split to get conformal splitting for the
        // prism. but for our application the tetrahedrons are only used for
        // inside detection so no need to get this information.
        //
        // I simply defined i2 as diagonal across from i6 and i3 diagonal
        // across from i5. and we can always split the pyramid as 
        // 1 (apex) 2 6 3 and 1 (apex) 2 6 5
        int i2=-1, i3=-1;
        for (int ii=0; ii<3; ii++) 
        {
            if (markedQuad[ii]) continue; 
            // at the non marked quad 
            for (int jj=0; jj<4; jj++) 
            {
                if (i6 == quadFaces(ii,jj)) 
                    i2 = quadFaces(ii,(jj+2)%4); 
                else if (i5 == quadFaces(ii,jj))
                    i3 = quadFaces(ii,(jj+2)%4); 
            }
        } 
        //std::cout << "Triangle set = \n" << triangleFaces << std::endl;
        //std::cout << "Quad set = \n" << quadFaces << std::endl;
        //std::cout << "tet A : " << i1 << ", " << i4 << ", " << i5 << ", " << i6 << std::endl;
        //std::cout << "tet B : " << i1 << ", " << i2 << ", " << i6 << ", " << i3 << std::endl;
        //std::cout << "tet C : " << i1 << ", " << i2 << ", " << i6 << ", " << i5 << std::endl;
        Tetrahedron tet_a(_pointPositions.row(i1), _pointPositions.row(i4), 
                          _pointPositions.row(i5), _pointPositions.row(i6)); 

        Tetrahedron tet_b(_pointPositions.row(i1), _pointPositions.row(i2), 
                          _pointPositions.row(i6), _pointPositions.row(i3)); 

        Tetrahedron tet_c(_pointPositions.row(i1), _pointPositions.row(i2), 
                          _pointPositions.row(i6), _pointPositions.row(i5)); 

        splittedTetrahedrons.push_back(tet_a); 
        splittedTetrahedrons.push_back(tet_b); 
        splittedTetrahedrons.push_back(tet_c); 
    }
    /* 
     *   top view: 
     *
     *  1     2
     *   o---o
     *   |\ /|
     *   | o | center point: 0
     *   |/ \|
     *   o---o
     *  4     3
     *
     *  without the need for face-orientation, the tetrahedrons can be 
     *  formed using (base face needs to be consistent) 
     *
     *  0 1 2 3, 0 1 4 3
     *
     */
    else if (N_faces==5 && N_overlappedPoints == 16) // pyramid case 
    {
        std::vector<int> pointIndicies(N_overlappedPoints); 
        Eigen::Vector4i baseFace; 
        int count = 0; 
        for (int ii=0; ii<N_faces; ii++)
        {
            const Eigen::VectorXi pointIndiciesThisFace = _faces[facesThisCellOwns[ii]]; 
            const int N_points = pointIndiciesThisFace.size(); 
            std::copy(pointIndiciesThisFace.data(), pointIndiciesThisFace.data()+N_points,pointIndicies.begin()+count);
            if (N_points==4) baseFace = pointIndiciesThisFace; 
            count += N_points; 
        }

        // decide which is apex by counting.. there should be a smarter way of
        // doing this
        std::vector<int> pointIndiciesTrimmed = STL_Wrapper::VectorSortAndTrim(pointIndicies); 
        int apexIndex=-1; 
        for (size_t ii=0; ii<pointIndiciesTrimmed.size(); ii++)
        {
            int indexOccurCount = std::count(pointIndicies.begin(), pointIndicies.end(), pointIndiciesTrimmed[ii]); 
            if (indexOccurCount == 4)
            {
                apexIndex = pointIndiciesTrimmed[ii]; 
                break; 
            } 
        } 
        STL_Wrapper::RemoveInPlace(pointIndiciesTrimmed, apexIndex); 

        // use apex index and base to construct two tetrahedrons
        Tetrahedron tet_a(_pointPositions.row(apexIndex ) , _pointPositions.row(baseFace(0)), 
                          _pointPositions.row(baseFace(1)), _pointPositions.row(baseFace(2))); 

        Tetrahedron tet_b(_pointPositions.row(apexIndex ) , _pointPositions.row(baseFace(0)), 
                          _pointPositions.row(baseFace(3)), _pointPositions.row(baseFace(2))); 

        splittedTetrahedrons.push_back(tet_a); 
        splittedTetrahedrons.push_back(tet_b); 
    }
    else
    {
        throw std::runtime_error("**ERROR** unknown non-tetrahedrons detected. don't know how to handle. face count: " + std::to_string(N_faces)+"; vertex count: "+std::to_string(N_overlappedPoints)); 
    }
}

void OpenfoamMesh::ReconstructTetrahedronMesh()
{
    std::cout << "reconstructing tetrahedron mesh" << std::endl; 

    const int N_cells = _cellIndicies.rows(); 
    const int N_owners = _owners.rows(); 
    const int N_neighbours = _neighbours.rows();

    std::vector<std::vector<int>> facesEachCellOwns(N_cells); 

    // loop through owners and neighbours to get the face count for each cell
    for (int ii=0; ii<N_owners; ii++) 
    {
        const int ind = _owners(ii) - _cellIndexOffset;
        if (ind>=0 && ind<N_cells)
            facesEachCellOwns[ind].push_back(ii); 
    }
    for (int ii=0; ii<N_neighbours; ii++) 
    {
        const int ind = _neighbours(ii) - _cellIndexOffset;
        if (ind>=0 && ind<N_cells)
            facesEachCellOwns[ind].push_back(ii); 
    }

    // create tetrahedrons and mark non-tetrahedron
    std::vector<int> nonTetrahedronCells;
    int N_tetrahedronsFromPyramids = 0; 
    int N_tetrahedronsFromPrisms   = 0; 
    for (int ii=0; ii<N_cells; ii++) 
    {
        // stored this information on the tet for later use
        TetrahedronExtraData extraData;
        extraData.openfoamCellIndex = ii + _cellIndexOffset;

        const int faceItHas = facesEachCellOwns[ii].size(); 
        if (faceItHas==4) 
        {
            std::vector<int> pointIndiciesForFace(faceItHas*3);  // if polygon has 4 faces, each face is triangle.
            for (int jj=0; jj<faceItHas; jj++) 
            {
                const int faceIndex   = facesEachCellOwns[ii][jj]; 
                const Eigen::VectorXi pointIndicies = _faces[faceIndex]; 
                std::copy(pointIndicies.data(), pointIndicies.data()+pointIndicies.size(), pointIndiciesForFace.begin()+jj*3);
            }
            // for our purpose no need for connectivity
            STL_Wrapper::VectorSortAndTrimInPlace(pointIndiciesForFace);
            Eigen::MatrixXd tetPointPositions(4,3); 
            for (int jj=0; jj<4; jj++)  // four vertex
                tetPointPositions.row(jj) = _pointPositions.row(pointIndiciesForFace[jj]); 

            Tetrahedron tet(tetPointPositions.row(0), tetPointPositions.row(1),
                            tetPointPositions.row(2), tetPointPositions.row(3)); 
            tet.SetExtraData(extraData); 
            AddTetrahedron(tet); 
        }
        else
        {
            nonTetrahedronCells.push_back(ii); 
            std::vector<Tetrahedron> splittedTetrahedrons; 
            SplitNonTetrahedrons(ii, facesEachCellOwns[ii], splittedTetrahedrons);
            if (splittedTetrahedrons.size()==3)
                N_tetrahedronsFromPrisms += 3; 
            else if (splittedTetrahedrons.size()==2)
                N_tetrahedronsFromPyramids += 2; 
            for (auto &tet : splittedTetrahedrons)
            {
                tet.SetExtraData(extraData);
                AddTetrahedron(tet); 
            }

        }
    }
    std::cout << " total of " << nonTetrahedronCells.size() << " cells are not tetrahedrons. \n"; 
    std::cout << "  number of tetrahedrons converted from prisms   : " << N_tetrahedronsFromPrisms   << std::endl;
    std::cout << "  number of tetrahedrons converted from pyramids : " << N_tetrahedronsFromPyramids << std::endl;
    std::cout << " total of " << N_tetrahedrons() << " tetrahedrons pushed into the tree mesh." << std::endl; 

    BuildTree(); 

    ///// visualize non-tetrahedrons points /////
    //Eigen::MatrixXd writeMatrix(nonTetrahedronCellPoints.size(), 3); 
    //for (size_t ii=0; ii<nonTetrahedronCellPoints.size(); ii++)
    //    writeMatrix.row(ii) = nonTetrahedronCellPoints[ii]; 
    //std::cout << writeMatrix << std::endl;
    //IO::writeMatrix_csv(writeMatrix,"data/nonTetrahedronCellPoints.csv"); 
}

void OpenfoamMesh::TrimAndCheckCellIndex()
{
    _cellIndexOffset = _cellIndicies.minCoeff(); 
    _cellIndicies.array() -= _cellIndexOffset; 

    for (int ii=0; ii<_cellIndicies.rows()-1; ii++) 
    {
        if (_cellIndicies(ii+1)-_cellIndicies(ii)!=1)
            throw std::runtime_error("**ERROR** cell ordering is wrong. cannot proceed"); 
    }

}

bool OpenfoamMesh::ReadFoamFile(const std::string &objectName, const int &N_cols, Eigen::MatrixXd &data)
{

    const std::string fileName = IO::AssembleFilePath(OpenfoamMesh::_polyMeshPath, objectName); 

    CheckValidity(fileName,1); 
    std::ifstream inFile(fileName); 

    const int N_readLines = ReadFoamHeader(inFile, objectName); 
    if (N_readLines>0)
        data = IO::readMatrixXd(N_readLines, N_cols, inFile, IO::BINARY);
    else 
        return false; 

    return true; 
}

bool OpenfoamMesh::ReadFoamFile_Int(const std::string &objectName, const int &N_cols, Eigen::MatrixXi &data, const std::string &zoneName) 
{

    const std::string fileName = IO::AssembleFilePath(OpenfoamMesh::_polyMeshPath, objectName); 

    CheckValidity(fileName,1); 
    std::ifstream inFile(fileName); 

    const int N_readLines = ReadFoamHeader(inFile, objectName, zoneName); 
    if (N_readLines>0)
        data = IO::readMatrixXi(N_readLines, N_cols, inFile, IO::BINARY);
    else 
        return false; 

    return true; 
}

bool OpenfoamMesh::ReadFoamFileFace(const std::string &objectName, CompressedList &data)
{
    const std::string fileName = IO::AssembleFilePath(OpenfoamMesh::_polyMeshPath, objectName); 

    CheckValidity(fileName,1); 
    std::ifstream inFile(fileName); 

    CheckFoamClass(inFile, "faceCompactList"); 

    Eigen::VectorXi &indiciesPointer = data.pointer; 
    Eigen::VectorXi &indicies = data.indicies; 

    const int N_compactFaceList = ReadFoamHeader(inFile, objectName); 

    // read compressed index
    indiciesPointer.resize(N_compactFaceList); 
    inFile.read((char*)indiciesPointer.data(), sizeof(int)*N_compactFaceList);

    // read the actual indicies after skipping the ascii stuff in between...
    const int N_indicies = indiciesPointer(N_compactFaceList-1);

    CheckReadBracket(inFile);  // ')'
    CheckReadFirstInteger(inFile, N_indicies);  // total number of indicies
    CheckReadBracket(inFile);  // '('

    indicies.resize(N_indicies);
    inFile.read((char*)indicies.data(), sizeof(int)*N_indicies); 

    if (indicies.size()==0 || indiciesPointer.size()==0)
        return false; 

    return true; 
}

///////////////////////////////////////////////////////////

void OpenfoamCase::ReinitializeCase()
{
    std::cout << "reinitializing openfoam case at path : " << _casePath << std::endl;
    SetProcessorPath(); 
    std::cout << " total of " << _processorPath.size() << " processor directories found" << std::endl;

    std::cout << " initializing underlying regular grid" << std::endl;
    Eigen::Vector3d minBound(-0.2178337591671500,-0.2017749603171500,-0.2192579190471500); 
    Eigen::Vector3d maxBound(0.2204221281671500,0.2364809270171500,0.2189979682871500); 
    Eigen::Vector3i cellCount(250,250,250); 

    _grid = UniformGrid<double>(minBound, maxBound, cellCount); 
}

void OpenfoamCase::SetProcessorPath()
{
    IO::listDirectoryMatch(_casePath.c_str(), "processor[[:digit:]]+", _processorPath);
    if (_processorPath.size()==0) throw std::runtime_error("**ERROR** could not find any processor directories in path : "+_casePath); 
    for (string &s : _processorPath) 
        s = IO::AssembleFilePath(_casePath, s); 
}


int OpenfoamCase::PrepareDataRead(const std::string &fieldName, const int &fieldCols, const double &dataStartTime, const double &dataStopTime) 
{

    _fieldCols = fieldCols; 
    _fieldName = fieldName; 

    std::string searchRangeString;
    if (dataStopTime<0) // not checking for stop time
        searchRangeString = "[" + std::to_string(dataStartTime) + ", inf)"; 
    else 
        searchRangeString = "[" + std::to_string(dataStartTime) + ", " + std::to_string(dataStopTime) + "]"; 

    std::cout << "searching all processor directories for data field : " << fieldName << " in the range " << searchRangeString << " ..." << std::endl;

    const bool checkAllProcessorDirectories = true; 

    if (N_processorDirectories()==0) 
        SetProcessorPath(); 

    const int N_processors = _processorPath.size(); 

    std::string &sampleDirectory = _processorPath[0];

    /* 
     * the data directories are always in the format
     * [0-9]*.[0-9]*
     */ 
    std::vector<std::string> timestepNames = IO::listDirectoryMatch(sampleDirectory.c_str(), "[[:digit:]]+\\.?[[:digit:]]*"); 

    // now check if field name is a file in all the time step directory
    for (std::string &t : timestepNames)
        CheckValidity(IO::AssembleFilePath(sampleDirectory,t), 1); 

    // check whether all the other processor directories have the same time step
    // folder and field. 
    if (checkAllProcessorDirectories) 
    {
        for (int ii=0; ii<N_processors; ii++) 
        {
            std::string &testDirectory = _processorPath[ii]; 
            std::vector<std::string> timestepNames_test = IO::listDirectoryMatch(testDirectory.c_str(), "[[:digit:]]+\\.?[[:digit:]]*");
            if (timestepNames!=timestepNames_test) throw std::runtime_error("**ERROR** time step range mismatched found in directory : "+testDirectory+" compared with the one in "+sampleDirectory);

            for (std::string &t : timestepNames)
                CheckValidity(IO::AssembleFilePath(testDirectory,t), 1); 
        } 
    }

    // remove everything not between data start time and data stop time
    std::vector<int> qualifiedTimestepsIndex; 
    for (size_t ii=0; ii<timestepNames.size(); ii++) 
    {
        const double timestep = std::atof(timestepNames[ii].c_str()); 
        if (timestep>=dataStartTime)
        {
            if (dataStopTime<0) // not checking for stop time
                qualifiedTimestepsIndex.push_back(ii); 
            else 
                if (timestep<=dataStopTime) 
                    qualifiedTimestepsIndex.push_back(ii); 
        }
    }

    // this will preserve the order
    _qualifiedTimesteps.reset(new std::vector<std::string>(qualifiedTimestepsIndex.size())); 
    _qualifiedTimesteps->resize(qualifiedTimestepsIndex.size());

    for (size_t ii=0; ii<qualifiedTimestepsIndex.size(); ii++)
        (*_qualifiedTimesteps)[ii] = timestepNames[qualifiedTimestepsIndex[ii]];


    // set processor data 
    for (int ii=0; ii<N_processors; ii++) 
    {
        std::shared_ptr<OpenfoamMesh> mesh(new OpenfoamMesh(_processorPath[ii], _zoneName)); 
        mesh->SetTimesteps(_qualifiedTimesteps); 
        mesh->SetFieldName(fieldName); 
        mesh->SetFieldCols(fieldCols);  

        _processorMesh.push_back(mesh); 
    }

    std::cout << " time steps ready to be read: \n "; 
    STL_Wrapper::PrintVectorContent(std::cout, *_qualifiedTimesteps, false, 20);

    return static_cast<int>(_qualifiedTimesteps->size()); 
}


void OpenfoamMesh::ReadAllTimesteps()
{
    _fieldData.resize(_timesteps->size()); 

    for (size_t ii=0; ii<_timesteps->size(); ii++)
    {
        std::string &str_timestep = (*_timesteps)[ii]; 
        std::ifstream inFile(IO::AssembleFilePath(_rootPath,str_timestep,_fieldName)); 
        int N_datapoints = ReadFoamHeader(inFile, _fieldName); 
        _fieldData[ii] = IO::readMatrixXd(N_datapoints, _fieldCols, inFile, IO::BINARY); 
    }
}

void OpenfoamMesh::GetCentroids(Eigen::MatrixXd &centroids)
{
    const int N_cells = _cellIndicies.rows(); 
    const int N_owners = _owners.rows(); 
    const int N_neighbours = _neighbours.rows();

    std::vector<std::vector<int>> facesEachCellOwns(N_cells); 

    centroids.resize(N_cells,3); 
    centroids.setZero();

    // loop through owners and neighbours to get the face count for each cell
    for (int ii=0; ii<N_owners; ii++) 
    {
        const int ind = _owners(ii) - _cellIndexOffset;
        if (ind>=0 && ind<N_cells)
            facesEachCellOwns[ind].push_back(ii); 
    }
    for (int ii=0; ii<N_neighbours; ii++) 
    {
        const int ind = _neighbours(ii) - _cellIndexOffset;
        if (ind>=0 && ind<N_cells)
            facesEachCellOwns[ind].push_back(ii); 
    }

    for (int ii=0; ii<N_cells; ii++) 
    {
        const int faceItHas = facesEachCellOwns[ii].size(); 

        std::vector<int> pointIndiciesForFace;  // if polygon has 4 faces, each face is triangle.
        for (int jj=0; jj<faceItHas; jj++)
        {
            const int faceIndex = facesEachCellOwns[ii][jj]; 
            const Eigen::VectorXi pointIndicies = _faces[faceIndex]; 
            for (int kk=0; kk<pointIndicies.size(); kk++) 
                pointIndiciesForFace.push_back(pointIndicies[kk]); 
        }
        STL_Wrapper::VectorSortAndTrimInPlace(pointIndiciesForFace);
        const int N_points = pointIndiciesForFace.size(); 
        Eigen::MatrixXd tetPointPositions(N_points,3); 
        for (int jj=0; jj<N_points; jj++)
            centroids.row(ii) += _pointPositions.row(pointIndiciesForFace[jj])/static_cast<double>(N_points);
    }

}

bool OpenfoamMesh::ReadTimestep(const std::string &timestep, const std::string &fieldName, const int &fieldCols, Eigen::MatrixXd &data) 
{
    //std::string &str_timestep = (*_timesteps)[ii]; 
    std::ifstream inFile(IO::AssembleFilePath(_rootPath, timestep, fieldName)); 
    if (!inFile) return false; 
    int N_datapoints = ReadFoamHeader(inFile, fieldName); 
    data = IO::readMatrixXd(N_datapoints, fieldCols, inFile, IO::BINARY); 
    return true; 
}



void OpenfoamCase::ReadAllTimesteps(const int &processorIndex, Eigen::MatrixXd &data) 
{

    std::cout << "reading data from all time steps" << std::endl;

    if (static_cast<int>(_processorMesh.size())<=processorIndex)
        throw std::runtime_error("**ERROR** processor index "+std::to_string(processorIndex)+" does not exist"); 

    for (auto & p : _processorMesh) 
        p->ReadAllTimesteps(); 
}



bool OpenfoamCase::ReadTimestep(const size_t &timeIndex, DataTimestep &data)
{
    if (timeIndex>=_qualifiedTimesteps->size()) return false;

    data.timestep  = _qualifiedTimesteps->at(timeIndex); 
    data.fieldName = _fieldName; 
    data.fieldCols = _fieldCols; 
    data.processorIndex.resize(_processorMesh.size()); 
    data.fieldData.resize(_processorMesh.size()); 

    for (size_t ii=0; ii<_processorMesh.size(); ii++)
    {
        data.processorIndex[ii] = ii; 
        bool success = _processorMesh[ii]->ReadTimestep(data.timestep,data.fieldName, data.fieldCols, data.fieldData[ii]); 
        if (!success)
            return false; 
    }

    return true;
}

void OpenfoamCase::GetMeshGridProjectionTable(const UniformGrid<double> &grid, const std::string &file_cachedIndicies, Eigen::MatrixXi &meshGridProjectionTable, const bool &readCached)
{

    std::cout << "constructing mesh-grid projection table. be patient..." << std::endl;
    if (meshGridProjectionTable.rows()!=0 || meshGridProjectionTable.cols()!=0)
    {
        std::cout << "**WARNING** projection table not empty, do nothing and keep the old one." << std::endl;
        return; 
    }

    const int N_processors = _processorMesh.size(); 

    // project it onto the grid 
    const Eigen::Vector3i cellCount = grid.GetCellCount(); 
    const int N_cells = cellCount[0] * cellCount[1] * cellCount[2]; 

    // this vector stores the tetrahedrons for each cell in the grid and the
    // processor directory it was found
    meshGridProjectionTable = -Eigen::MatrixXi::Ones(N_cells,2);  // cellindicies, processorindicies
    std::string file_meshGridProjectionTable(file_cachedIndicies); 

    int count =0;

    if (IO::ExistFile(file_meshGridProjectionTable) && readCached)
    {
        std::cout << "  read from hard drive : " << file_meshGridProjectionTable << std::endl;
        IO::readMatrixX<int>(meshGridProjectionTable, file_meshGridProjectionTable.c_str(), IO::BINARY);
    }
    else
    {
        #pragma omp parallel for
        for (int kk=0; kk<cellCount[2]; kk++) 
        {
            for (int jj=0; jj<cellCount[1]; jj++) 
            {
                for (int ii=0; ii<cellCount[0]; ii++) 
                {
                    const Eigen::Vector3d cellPosition = grid.GetCellCenterPosition(ii,jj,kk); 
                    const int index = FlattenIndex(cellCount[0],cellCount[1],cellCount[2],ii,jj,kk); 
                    if (cellPosition.norm()>0.2) continue; // outside the spherical region 

                    for (int pp=0; pp<N_processors; pp++) 
                    {
                        std::shared_ptr<OpenfoamMesh> mesh = _processorMesh[pp]; 
                        const int tetIndex = mesh->FindEnclosedTetrahedron(cellPosition,50); 
                        if (tetIndex >0) 
                        {
                            meshGridProjectionTable(index,0) = mesh->Get(tetIndex).GetExtraData().openfoamCellIndex; 
                            meshGridProjectionTable(index,1) = pp; 
                            break;
                        }
                    }
                    #pragma omp critical 
                    {
                        count ++; 
                        if (count %(10000) ==0)
                            std::cout << "  computed " << count << " grid points. \r" << std::flush;
                    } // pragma omp critical
                }
            }
        }
        std::cout << std::endl;
        IO::writeMatrixX<int>(meshGridProjectionTable, file_meshGridProjectionTable.c_str(), IO::BINARY); 
    }
}

void OpenfoamCase::ProjectDataTimestep(const DataTimestep &data, const Eigen::MatrixXi &meshGridProjectionTable, const UniformGrid<double> &grid, Eigen::MatrixXd &projectedData)
{
    const Eigen::Vector3i cellCount = grid.GetCellCount(); 
    const int N_cells = cellCount[0]*cellCount[1]*cellCount[2];
    projectedData.setZero(N_cells, data.fieldCols); 

    for (int kk=0; kk<cellCount[2]; kk++) 
    {
        for (int jj=0; jj<cellCount[1]; jj++)
        {
            for (int ii=0; ii<cellCount[0]; ii++)
            {
                const int index  = FlattenIndex(cellCount[0],cellCount[1],cellCount[2],ii,jj,kk); 
                const int cellIndex = meshGridProjectionTable(index,0); 
                const int processorIndex = meshGridProjectionTable(index,1); 

                if (cellIndex<0 || processorIndex<0) // enclosed tet not found 
                    continue; 

                projectedData.row(index) = data.fieldData[processorIndex].row(cellIndex); 
            }
        }
    }
}

// this data has to be in the flatten format prepared by projectDataTimestep
void OpenfoamCase::WriteDataTimestepVTK(const std::string &fileName, const std::string &fieldName, const Eigen::MatrixXd &data, const UniformGrid<double> &grid)
{
    const Eigen::Vector3i cellCount = grid.GetCellCount(); 
    Eigen::MatrixXd gridPosition(cellCount[0]*cellCount[1]*cellCount[2],3); 

    for (int kk=0; kk<cellCount[2]; kk++) 
    {
        for (int jj=0; jj<cellCount[1]; jj++)
        {
            for (int ii=0; ii<cellCount[0]; ii++)
            {
                const int index  = FlattenIndex(cellCount[0],cellCount[1],cellCount[2],ii,jj,kk); 
                gridPosition.row(index) = _grid.GetCellCenterPosition(ii,jj,kk);
            }
        }
    }

    VTKConverter::VTKStructureGridWithScalarFromEigen( gridPosition, data, fileName, "fieldName", VTKConverter::BINARY, cellCount );

}

void OpenfoamCase::GetCentroids(Eigen::MatrixXd &centroids) 
{

    //centroids = Eigen::MatrixXd(); 
    const int N_processors = _processorMesh.size(); 
    for (int ii=0; ii<N_processors; ii++) 
    {
        Eigen::MatrixXd centroidsProc; 
        _processorMesh[ii]->GetCentroids(centroidsProc);
        const int oldSize = centroids.rows(); 
        centroids.conservativeResize(centroidsProc.rows()+oldSize,3);
        centroids.block(oldSize,0,centroidsProc.rows(),3) = centroidsProc; 
    }

}

