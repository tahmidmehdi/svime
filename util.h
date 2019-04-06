/*
 * util.h
 * Utilities
 *  Updated on: Feb 21, 2019
 *      Author: Tahmid Mehdi
 *
Copyright 2019 Tahmid Mehdi
This file is part of svime.

svime is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

svime is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with svime.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace Util {
// Save matrix to file
template<class Matrix>
void write_binary(const char* filename, const Matrix& matrix){
    std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(typename Matrix::Index));
    out.write((char*) (&cols), sizeof(typename Matrix::Index));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
    out.close();
}
// Load matrix to file
template<class Matrix>
void read_binary(const char* filename, Matrix& matrix){
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    typename Matrix::Index rows=0, cols=0;
    in.read((char*) (&rows),sizeof(typename Matrix::Index));
    in.read((char*) (&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read((char*) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
    in.close();
}
// Get index of max element in each row of matrix
template<class Matrix>
std::vector<typename Matrix::Index> max_index(Matrix& matrix){
	int n = matrix.rows();
	std::vector<typename Matrix::Index> maxIndex(n, 0);
	std::vector<float> maxVal(n, 0);
	for( int i = 0; i < n; ++i ) {
	    maxVal[i] = matrix.row(i).maxCoeff( &maxIndex[i] );
	}
	return maxIndex;
}
// Load sequences from filename to seqs
void loadSeqs(std::string filename, std::vector<std::string>& seqs) {
	seqs.clear();
	std::ifstream file (filename.c_str());
	std::string line;
	while ( std::getline (file, line) )
	{
	   seqs.push_back(line);
	}
}
// Access v[start:end]
template<typename T>
std::vector<T> slice(std::vector<T> const& v, int start, int end) {
	auto first = v.cbegin()+start;
	auto last = v.cbegin()+end;
	std::vector<T> vec(first, last);
	return vec;
}

} // Util::


#endif /* SRC_UTIL_H_ */
