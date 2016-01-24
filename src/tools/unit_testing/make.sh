rm OpenfoamMesh.o TetrahedronMesh.o test
g++ -Wall -O3 -std=c++11 -c TetrahedronMesh.cpp -I../compute-quadrupole/include -I/home/jui-hsien/code/turbsound_postfluid/tools/IO
g++ -Wall -O3 -std=c++11 -c OpenfoamMesh.cpp -I../compute-quadrupole/include -I/home/jui-hsien/code/turbsound_postfluid/tools/IO -DUSE_BOOST -lboost_regex -fopenmp
g++ -Wall -O3 -std=c++11 test.cpp -o test OpenFOAM_Helper.o OpenfoamMesh.o TetrahedronMesh.o -I../compute-quadrupole/include -I../../tools/IO -lboost_timer -lboost_system -lboost_regex -DUSE_BOOST -fopenmp
