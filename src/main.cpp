/*
 * main.cpp
 *
 *  Created on: Jul 13, 2018
 *      Author: Sarah Lutteropp, Niklas Birth
 */

#include <iostream>
#include <string>
#include <array>
#include <fstream>

#include "quartet_topology_checker.hpp"

// Adapted from https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/14266139#14266139
// Why does C++ not have anything like this included in its standard library? 
std::vector<std::string> split_string(std::string s, std::string delimiter) {
	size_t pos = 0;
	std::vector<std::string> v;
	std::string token;

	while ((pos = s.find(delimiter)) != std::string::npos) {
	    v.push_back(s.substr(0, pos));
	    s.erase(0, pos + delimiter.length());
	}

	v.push_back(s);

	return v;
}

void evaluateQuartets(const std::string& fastaPath, const std::string& speciesTreePath, const std::string& quartetsPath) {
	TopologyChecker topoChecker(fastaPath, speciesTreePath);

	std::ifstream infile(quartetsPath);
	std::string qtree;
	while (std::getline(infile, qtree)) {
		qtree = qtree.substr(2, qtree.length() - 5);
		std::vector<std::string> v = split_string(qtree, ",");

		int a = topoChecker.labelToInt(v[0]);
		int b = topoChecker.labelToInt(v[1].substr(0, v[1].length() - 1));
		int c = topoChecker.labelToInt(v[2].substr(1, v[2].length()));
		int d = topoChecker.labelToInt(v[3]);

		std::cout << topoChecker.sameTopologyAsReference(a, b, c, d) << "\n";
	}
}

int main(int argc, char* argv[]) {
	if (argc != 4) {
		std::cout << "Usage: ./single_quartet_check <FASTA file> <quartet tree file> <reference tree file>\n";
		return 1;
	}

	std::string fastaPath = argv[1];
	std::string quartetsPath = argv[2];
	std::string treePath = argv[3];

	evaluateQuartets(fastaPath, treePath, quartetsPath);
}