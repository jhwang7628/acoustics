find_path(Eigen_INCLUDE_DIR Eigen/Core                                                                                      
  HINTS "${EIGEN_ROOT}"                                                                                                
  HINTS "/usr/local/Cellar/"
  # HINTS "/home/jw969/opt/eigen"
  PATH_SUFFIXES "eigen"                                                                                                
  DOC "The directory where Eigen/Core resides"                                                                                                
)                                                                                                                                             
                                                                                                                                              
mark_as_advanced(Eigen_INCLUDE_DIR)                                                                                                           
                                                                                                                                              
include(FindPackageHandleStandardArgs)                                                                                                        
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Eigen DEFAULT_MSG Eigen_INCLUDE_DIR) 
