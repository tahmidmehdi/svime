/*
 * processFasta.cpp
 *
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

#include "processFasta.h"
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <boost/algorithm/string.hpp>
#include <iostream>

void saveData(std::string filename, const std::vector<std::string>& data) {
	std::ofstream outFile(filename);
	for (const auto &e : data) outFile << e << "\n";
}

std::map<std::string, int> faToMatrix(std::string fa, int w, std::string outDir) {
	std::vector<std::string> regions;
	std::vector<std::string> seqs;
	std::string line;
	std::ifstream faFile (fa.c_str());
	std::string delimiters = ":-";
	std::vector<std::string> parts;
	std::string chr, seqFile, regionsFile;
	std::map<std::string, int> regionsPerChr;
	int start, end;
	std::string prevChr = "chr1";

	if (faFile.is_open()) {
		while ( std::getline (faFile, line) ) {
			if (line[0] == '>') { // coordinate line
				line.erase(line.begin());
				// split line to extract chromosome, start & end coordinates
				boost::split(parts, line, boost::is_any_of(delimiters));
				chr = parts[0]; start = std::stoi(parts[1]); end = std::stoi(parts[2]);
				if (chr != prevChr) { // if new chromosome
					seqFile = outDir+"/"+prevChr+"_seq.txt";
					regionsFile = outDir+"/"+prevChr+"_coords.txt";
					saveData(seqFile, seqs);
					saveData(regionsFile, regions);
					// store chromosome & number of w-mers in dictionary
					regionsPerChr.insert(std::pair<std::string, int>(prevChr, regions.size()));
					seqs.clear(); regions.clear();
				}
				prevChr = chr;
			} else {
				boost::to_upper(line);
				// store sliding windows
				for( int i = 0; i <= end-start-w; ++i ) {
					seqs.push_back (line.substr(i, w));
					regions.push_back (chr+":"+std::to_string(start+i)+"-"+std::to_string(start+i+w));
				}
			}
	    }
		faFile.close();
	}
	return regionsPerChr;
}
