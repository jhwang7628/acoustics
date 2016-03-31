#include <linearalgebra/SparseLinearSystemSolver.h> 
#include <utils/STL_Wrapper.h> 

#include <Eigen/IterativeLinearSolvers>


void SparseLinearSystemSolver::StageEntry(const int &x_i, const int &x_j, const double &e_ij)
{
    if (x_i >= _N || x_j >= _N) throw std::runtime_error("**ERROR** requested node out-of-bounds"); 
    _stagedEntries.push_back(Triplet(x_i, x_j, e_ij)); 

    _stagedEntriesPattern.push_back(TripletInt(x_i, x_j, 1)); 
    _stagedEntriesPattern.push_back(TripletInt(x_j, x_i, 1)); 
}

void SparseLinearSystemSolver::FillIn(const bool &destroyOld)
{
    if (!destroyOld) 
    {
        // fill-in for A
        for (int row=0; row<_A.outerSize(); row++) 
            for (SparseMatrix::InnerIterator it(_A,row); it; ++it) 
                _stagedEntries.push_back(Triplet(it.row(), it.col(), it.value())); 

        // fill-in for A patter
        for (int row=0; row<_A_pattern.outerSize(); row++) 
            for (SparseMatrixPattern::InnerIterator it(_A_pattern,row); it; ++it) 
                _stagedEntriesPattern.push_back(TripletInt(it.row(), it.col(), it.value())); 
    }

    const int addedEntries = _stagedEntries.size();
    if (addedEntries==0) return; 

    _A.setFromTriplets(_stagedEntries.begin(), _stagedEntries.end()); 
    _stagedEntries.clear(); 

    _A_pattern.setFromTriplets(_stagedEntriesPattern.begin(), _stagedEntriesPattern.end()); 

    _M = addedEntries; 

}

void SparseLinearSystemSolver::GetAllNonzeroElements(Eigen::MatrixXd &nonzeroElements)
{

    nonzeroElements.resize(_A.nonZeros(), 3); 

    int count = 0;
    for (int row=0; row<_A.outerSize(); row++) 
        for (SparseMatrix::InnerIterator it(_A,row); it; ++it) 
        {
            nonzeroElements(count,0) = it.row(); 
            nonzeroElements(count,1) = it.col(); 
            nonzeroElements(count,2) = it.value(); 

            count ++; 
        } 
}

void SparseLinearSystemSolver::GetAllComponents(std::vector<std::vector<int> > &allComponents)
{
    _belongComponent.clear(); 
    _belongComponent.resize(_N,-1); 

    allComponents.clear(); 

    int root; 

    for (int row=0; row<_A_pattern.outerSize(); row++) 
    {
        if (_A_pattern.row(row).nonZeros() != 0)
        {
            // use the first element in the column as the root for BFS 
            SparseMatrixPattern::InnerIterator it(_A_pattern,row); 
            root = it.col(); 

            if (_belongComponent[root] == -1)  // not searched yet
            {
                std::vector<int> component; 
                GetComponent(root, component, allComponents.size()); 
                allComponents.push_back(component); 
            } 
        }
    }
}


// BFS graph traversal 
// Ref: https://en.wikipedia.org/wiki/Breadth-first_search
void SparseLinearSystemSolver::GetComponent(const int &root, std::vector<int> &component, const int &componentID)
{

    if (root >= _N) throw std::runtime_error("**ERROR** given root is larger than number of nodes in the graph"); 

    component.clear();

    std::vector<int> distances(_N, std::numeric_limits<int>::max()); 
    std::vector<int>   parents(_N, -1); 

    distances[root] = 0;
    component.push_back(root); 
    _belongComponent[root] = componentID; 
    std::queue<int> searchQueue; 
    searchQueue.push(root); 

    int current; 
    int currentNeighbour; 

    while(!searchQueue.empty())
    {
        current = searchQueue.front(); 
        searchQueue.pop(); 
        for (SparseMatrixPattern::InnerIterator it(_A_pattern,current); it; ++it)
        {
            currentNeighbour = it.col(); 
            if (distances[currentNeighbour] == std::numeric_limits<int>::max()) // not visited yet
            { 

                _belongComponent[currentNeighbour] = componentID; 

                component.push_back(currentNeighbour); 
                distances[currentNeighbour] = distances[current] + 1; 
                parents[currentNeighbour] = current; 
                searchQueue.push(currentNeighbour); 
            }
        }
        
    }

}

void SparseLinearSystemSolver::DecoupledSolve(Eigen::VectorXd &x) 
{

    std::vector<std::vector<int> > allComponents; 
    GetAllComponents(allComponents); 

    const int numberComponents = allComponents.size(); 

    Eigen::MatrixXd sub_A; 
    Eigen::VectorXd sub_b; 
    Eigen::VectorXd sub_x; 

    x.resize(_N); 

    for (int cc=0; cc<numberComponents; cc++) 
    {

        std::vector<int> &component = allComponents[cc]; 
        std::sort(component.begin(), component.end()); 
        const int subProblemSize = component.size(); 
        sub_A.setZero(subProblemSize, subProblemSize); 
        sub_b.setZero(subProblemSize); 

        {
            //boost::timer::auto_cpu_timer t("old method %w sec\n"); 
            // fill the linear system 
            // TODO this might be inefficient due to bisection random access
            for (int ii=0; ii<subProblemSize; ii++) 
                for (int jj=0; jj<subProblemSize; jj++) 
                    sub_A(ii,jj) = _A.coeff(component[ii],component[jj]);
        } 

        //std::cout << " component = \n";
        //STL_Wrapper::PrintVectorContent(std::cout, component); 

        //std::cout << "(old) sub_A = \n" << sub_A << std::endl;
        //sub_A.setZero();

        //{
        //    boost::timer::auto_cpu_timer t("new method %w sec\n"); 




        //    //for (int ii=0; ii<subProblemSize; ii++) 
        //    //{
        //    //    SparseMatrix::InnerIterator it(_A,component[ii]); 

        //    //    for (int jj=0; jj<subProblemSize; jj++) 
        //    //    {
        //    //        if (component[jj]==it.col()) // match
        //    //        {
        //    //            sub_A(ii,jj) = it.value(); 
        //    //            ++it; 
        //    //        }
        //    //        else 
        //    //        {
        //    //            sub_A(ii,jj) = 0; 
        //    //        }
        //    //    }

        //    //}
        //} 

        //std::cout << "(new) sub_A = \n" << sub_A << std::endl;

        for (int ii=0; ii<subProblemSize; ii++) 
            sub_b(ii) = _bDense(component[ii]); 

        sub_x = sub_A.householderQr().solve(sub_b); 

        // fill back x 
        for (int ii=0; ii<subProblemSize; ii++) 
            x(component[ii]) = sub_x(ii);

    }
}

void SparseLinearSystemSolver::SparseSolve(Eigen::VectorXd &x)
{
    _A.makeCompressed();

    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> solver; 
    solver.compute(_A); 

    //Eigen::SparseLU<SparseMatrix> solver; 
    //solver.analyzePattern(_A);
    //solver.factorize(_A); 


    //Eigen::VectorXd bDense = Eigen::VectorXd(_b); 

    //std::cout << "bDense = " << bDense;

    x = solver.solve(_bDense); 
}

void SparseLinearSystemSolver::BiCGStabSolve(Eigen::VectorXd &x) 
{
    Eigen::BiCGSTAB<SparseMatrix> solver; 
    solver.compute(_A); 

    //Eigen::VectorXd bDense = Eigen::VectorXd(_b); 

    x = solver.solve(_bDense); 
    
}

void SparseLinearSystemSolver::DenseSolve(Eigen::VectorXd &x) 
{
    Eigen::MatrixXd denseMatrix = Eigen::MatrixXd(_A); 
    //Eigen::VectorXd denseRHS    = Eigen::VectorXd(_b); 

    x = denseMatrix.householderQr().solve(_bDense);
}

void SparseLinearSystemSolver::PrintMatrixDense(std::ostream &os)
{
    Eigen::MatrixXd denseMatrix = Eigen::MatrixXd(_A); 
    os << "A=\n" << denseMatrix << std::endl;
}

void SparseLinearSystemSolver::PrintVectorDense(std::ostream &os)
{
    Eigen::VectorXd denseVector = Eigen::VectorXd(_b); 
    os << "b=\n" << denseVector << std::endl;
}


void SparseLinearSystemSolver::GetAllComponentsVerbose(std::vector<std::vector<int> > &allComponents)
{
    GetAllComponents(allComponents); 

    std::cout << "there are " << allComponents.size() << " connected components found in the graph" << std::endl; 
    for (size_t ii=0; ii<allComponents.size(); ii++) 
    {
        std::cout << "component " << ii << " : "; 
        STL_Wrapper::PrintVectorContent(std::cout, allComponents[ii]); 
    }



}

