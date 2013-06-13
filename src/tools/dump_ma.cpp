/*
 * =====================================================================================
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  dump_ma.cpp
 *
 *        Version:  1.0
 *        Created:  12/12/11 16:43:50
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include <iostream>
#include "io/MatrixIO.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    if ( argc != 2 ) 
    {
        cerr << "Invalid arguments!" << endl;
        return 1;
    }

    Matrix<double>* mat = load_ma_matrixd(argv[1]);
    const Matrix<double>::size_type* s = mat->shape();
    cout << "SIZE: " << s[0] << ' ' << s[1] << endl;
    for(int i = 0;i < s[0];++ i)
    {
        for(int j = 0;j < s[1];++ j)
            cout << (*mat)[i][j] << ' ';
        cout << endl;
    }
    return 0;
}
