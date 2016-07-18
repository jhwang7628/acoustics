#include <Eigen/Dense> 
#include <iostream>

int main()
{
    Eigen::VectorXd vec(5); 
    vec << 1, 2, 3, 4, 5; 

    std::cout << vec << std::endl;
    std::cout << &vec << std::endl;

    // this does not work, eigen will create copy
    Eigen::VectorXd vec_noCopy; 
    vec_noCopy = *(&vec); 

    std::cout << vec << std::endl;
    std::cout << &vec << std::endl;

    return 0;
}
