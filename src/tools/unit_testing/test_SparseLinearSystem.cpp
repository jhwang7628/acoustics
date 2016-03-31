#include <iostream> 
#include <fstream>

#include <utils/STL_Wrapper.h>
#include <tools/unit_testing/testing.h> 
#include <linearalgebra/SparseLinearSystemSolver.h> 


#include <boost/timer/timer.hpp> 


int main(int argc, char **argv) {


    if (argc < 3) 
    {
        std::cerr << "**Usage: " << argv[0] << " <problem_dimension> <number_nonzero_element> [test_with_dense_solve] [output_text_file]\n"; 
        exit(1); 
    }


    bool testDense; 

    const int dim = atoi(argv[1]); 
    const int nnz = atoi(argv[2]);
    char *outFile; 
    if (argc >= 5)
        outFile = argv[4]; 
    else 
        outFile = nullptr; 


    if (argc >= 4) 
        testDense = (atoi(argv[3])==1 ? true : false); 
    else 
        testDense = true; 

    std::cout << "testing on SparseLinearSystemSolverclass for dim = " << dim << "; nnz = " << nnz << " \n"; 

    //const int dim = 10000; 
    //const int nnz = 500;
    SparseLinearSystemSolver spLinear(dim); 



    srand(time(NULL)); 
    for (int ii=0; ii<nnz; ii++) 
    {

        const int randomRow = rand() % dim; 
        const int randomCol = rand() % dim; 

        const double randomEntry = (double)rand() / (double)RAND_MAX; 

        spLinear.StageEntry(randomRow, randomRow, 1.0); 
        spLinear.StageEntry(randomCol, randomCol, 1.0);
        spLinear.StageEntry(randomRow, randomCol, randomEntry); 


        const double randomRHS1 = (double)rand() / (double)RAND_MAX; 
        const double randomRHS2 = (double)rand() / (double)RAND_MAX; 
        spLinear.FillInRHS(randomRow, randomRHS1); 
        spLinear.FillInRHS(randomCol, randomRHS2); 

        
    }

    for (int ii=0; ii<dim; ii++) spLinear.StageEntry(ii,ii,1.0); 

    //spLinear.StageEntry(0,1,2); 
    //spLinear.StageEntry(1,0,3); 
    //spLinear.StageEntry(4,3,4); 

    //spLinear.FillInRHS(1,0.5); 
    //spLinear.FillInRHS(2,1.0); 
    //spLinear.FillInRHS(3,2.0); 
    //spLinear.FillInRHS(8,2.0); 

    //spLinear.StageEntry(0,1,1.0); 
    //spLinear.StageEntry(0,2,1.0); 
    //spLinear.StageEntry(1,2,1.0); 

    //spLinear.StageEntry(3,4,1.0); 
    //spLinear.StageEntry(4,5,1.0); 
    //spLinear.StageEntry(5,6,1.0); 
    //spLinear.StageEntry(4,6,1.0); 

    //spLinear.StageEntry(8,9,1.0);
    //spLinear.StageEntry(8,1,1.0);



    //spLinear.AddEdgeToStage(40,5); 
    //spLinear.StageEntry(39,1,1.0); 

    spLinear.FillIn(); 

    Eigen::MatrixXd entries; 
    spLinear.GetAllNonzeroElements(entries); 

    //COUT_SDUMP(entries); 

    //spLinear.PrintMatrixDense(std::cout);
    //spLinear.PrintVectorDense(std::cout); 

    Eigen::VectorXd x; 
    Eigen::VectorXd xDense; 
    Eigen::VectorXd xBiCGStab; 
    Eigen::VectorXd xDecoupled; 
    Eigen::VectorXd xExact; 

    std::ostream *os; 
    if (outFile) 
        os = new std::ofstream(outFile, std::fstream::app); 
    else 
        os = &std::cout; 

    if ( (int)os->tellp() == 0 ) 
    {
        if (testDense) 
            *os << "dimension nnz dense_solve_time sparse_solve_time BiCGStab_solve_time decoupled_solve_time sparse_error BiCGStab_error decoupled_error\n"; 
        else 
            *os << "dimension nnz sparse_solve_time BiCGStab_solve_time decoupled_solve_time sparse_error BiCGStab_error decoupled_error\n"; 
    }


    *os << dim << " " << nnz << " "; 


    const bool test_sparseQR = false; 
    const bool test_BiCGStab = true; 
    const bool test_decoupled = false; 

    if (testDense)
    {
        boost::timer::auto_cpu_timer t(*os, 16, "%w "); 
        spLinear.DenseSolve(xDense); 
    }

    if (test_sparseQR)
    {
        boost::timer::auto_cpu_timer t(*os, 16, "%w "); 
        spLinear.SparseSolve(x); 
    }

    if (test_BiCGStab)
    {
        boost::timer::auto_cpu_timer t(*os, 16, "%w "); 
        spLinear.BiCGStabSolve(xBiCGStab); 
    }

    if (test_decoupled)
    {
        boost::timer::auto_cpu_timer t(*os, 16, "%w "); 
        spLinear.DecoupledSolve(xDecoupled); 
    }


    if (testDense) 
        xExact = xDense; 
    else 
        if (test_sparseQR)
            xExact = x; 
        else 
            xExact = xBiCGStab; 




    //std::cout << "sparse    solve error = " << (x          - xDense).norm()/xDense.norm() << std::endl; 
    //std::cout << "decoupled solve error = " << (xDecoupled - xDense).norm()/xDense.norm() << std::endl; 

    if (test_sparseQR)
        *os << (x          - xExact).norm()/xExact.norm() << " "; 
    if (test_BiCGStab)
        *os << (xBiCGStab  - xExact).norm()/xExact.norm() << " "; 
    if (test_decoupled)
        *os << (xDecoupled - xExact).norm()/xExact.norm() << " "; 

    *os << std::endl; 
    if (outFile)
        reinterpret_cast<std::ofstream*>(os)->close(); 


    return 0; 
}
