#ifndef OPENFOAMMESH_H 
#define OPENFOAMMESH_H 

#include <Eigen/Dense> 

#include "utils/IO/IO.h" 
#include "TetrahedronMesh.h" 
#include "Grid.h" 


/* 
 * similar to compressed sparse matrix storage CRS
 */
struct CompressedList
{
    Eigen::VectorXi     pointer; 
    Eigen::VectorXi     indicies; 

    inline size_t Size() const { return pointer.size()-1; } 
    inline size_t N_polygons(const size_t ind) const { return pointer(ind+1)-pointer(ind); }

    Eigen::VectorXi operator[] (const size_t ind)
    {
        if (ind>=Size()) throw std::runtime_error("**ERROR** index out of bounce : "+std::to_string(ind)+"->"+std::to_string(Size())); 
        return indicies.segment(pointer(ind),pointer(ind+1)-pointer(ind)); 
    }
};

struct DataTimestep
{
    std::string                     timestep; 
    std::string                     fieldName; 
    int                             fieldCols; 

    std::vector<int>                processorIndex; 
    std::vector<Eigen::MatrixXd>    fieldData;  // for each processor directories (sorted by the accessed order)
};


class OpenfoamMesh : public TetrahedronMesh
{
    private: 
        std::string         _rootPath; 
        std::string         _polyMeshPath; 
        std::string         _zoneName; 

        Eigen::MatrixXd     _pointPositions; 
        Eigen::MatrixXi     _cellIndicies; 
        Eigen::MatrixXi     _owners; 
        Eigen::MatrixXi     _neighbours; 
        CompressedList      _faces; 

        int                 _cellIndexOffset; // the zone we are looking at might start from nonzero index

        int                                         _fieldCols; 
        std::string                                 _fieldName; 
        std::vector<Eigen::MatrixXd>                _fieldData;  // matrixxd for each time step 
        std::shared_ptr<std::vector<std::string>>   _timesteps; 
        

    public: 

        OpenfoamMesh(){ }
        OpenfoamMesh(const std::string &root, const std::string &zone) : 
            _rootPath(root),
            _polyMeshPath(IO::AssembleFilePath(root,"constant/polyMesh")), 
            _zoneName(zone){
                ReinitializeMesh(); 
                ReconstructTetrahedronMesh(); 
            }

        inline void SetTimesteps(std::shared_ptr<std::vector<std::string>> p) {_timesteps = p;} 
        inline void SetFieldName(const string &fieldName) {_fieldName=fieldName;}
        inline void SetFieldCols(const int &fieldCols) {_fieldCols=fieldCols;} 
        inline int GetCellIndexOffset(){ return _cellIndexOffset; }

        inline Eigen::MatrixXd& GetFieldData(const int &timeIndex){return _fieldData[timeIndex];}
        //inline std::string                               GetFieldName(){return _fieldName;} 
        //inline std::string                               GetFieldCols(){return _fieldCols;} 
        //inline std::shared_ptr<std::vector<std::string>> GetTimesteps()(return _timesteps;} 

        void ReinitializeMesh();
        int CountNonTriangleFaces() const; 
        void SplitNonTetrahedrons(const int &nonTetrahedronCellIndex, const std::vector<int> &facesThisCellOwns, std::vector<Tetrahedron> &); 
        void ReconstructTetrahedronMesh();
        void TrimAndCheckCellIndex(); 
        bool ReadFoamFile    (const std::string &, const int &, Eigen::MatrixXd &); 
        bool ReadFoamFile_Int(const std::string &, const int &, Eigen::MatrixXi &, const std::string &zoneName="NO_ZONE");
        bool ReadFoamFileFace(const std::string &, CompressedList &);

        bool ReadTimestep(const std::string &timestep, const std::string &fieldName, const int &fieldCols, Eigen::MatrixXd &data);
        void ReadAllTimesteps(); 
        void GetCentroids(Eigen::MatrixXd &centroids); 

};

class OpenfoamCase
{
    public: 

    private: 
        std::string                                 _casePath; 
        std::string                                 _zoneName; 
        std::vector<std::string>                    _processorPath; 
        std::shared_ptr<std::vector<std::string>>   _qualifiedTimesteps; 
        std::vector<std::shared_ptr<OpenfoamMesh>>  _processorMesh; 

        int                                         _fieldCols; 
        std::string                                 _fieldName; 

        UniformGrid                                 _grid; 


    public: 
        inline int N_processorDirectories() { return _processorPath.size(); } 

        void ReinitializeCase(); 
        void SetProcessorPath();
        int PrepareDataRead(const std::string &fieldName, const int &fieldCols, const double &dataStartTime=0, const double &dataStopTime=-1); 
        void ReadAllTimesteps(const int &processorIndex, Eigen::MatrixXd &data); 
        bool ReadTimestep(const size_t &timeIndex, DataTimestep &data);
        void GetMeshGridProjectionTable(const UniformGrid &grid,const std::string &file_cachedIndicies,Eigen::MatrixXi &meshGridProjectionTable, const bool &readCached); 
        void ProjectDataTimestep(const DataTimestep &data, const Eigen::MatrixXi &meshGridProjectionTable, const UniformGrid &grid, Eigen::MatrixXd &projectedData); 
        void WriteDataTimestepVTK(const std::string &fileName, const std::string &fieldName, const Eigen::MatrixXd &data, const UniformGrid &grid); 

        inline int GetCellIndexOffsetProcessor(const int &proc){ return _processorMesh.at(proc)->GetCellIndexOffset(); }
        void GetCentroids(Eigen::MatrixXd &centroids); 

        OpenfoamCase(const std::string &casePath, const std::string &zone) : 
            _casePath(casePath),
            _zoneName(zone) {
            }

}; 

#endif
