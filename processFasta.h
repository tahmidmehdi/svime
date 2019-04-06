/*
 * processFasta.h
 * Process FASTA files
 *  Updated on: Feb 16, 2019
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

#ifndef SRC_PROCESSFASTA_H_
#define SRC_PROCESSFASTA_H_
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Save sequences or genomic coordinates contained in data to a text file named filename
void saveData(std::string filename, const std::vector<std::string>& data);
/* Saves 2 files for each chromosome in fa (a FASTA file): One file that stores
all overlapping w-mers and another that stores coordinates for each w-mer. It
produces a map where each key is a chromosome and its value is the number of
w-mers in the chromosome. */
std::map<std::string, int> faToMatrix(std::string fa, int w, std::string outDir);

#endif /* SRC_PROCESSFASTA_H_ */
