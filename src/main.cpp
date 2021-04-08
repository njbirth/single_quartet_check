/*
 * main.cpp
 *
 *  Created on: Jul 13, 2018
 *      Author: Sarah Lutteropp
 */

#include <iostream>
#include <string>
#include <array>

#include "quartet_topology_checker.hpp"

void evaluateQuartets(const std::string& fastaPath, const std::string& speciesTreePath, const std::string& multiSPAMPath) {
	size_t correct = 0;
	size_t wrong = 0;

	TopologyChecker topoChecker(fastaPath, speciesTreePath, multiSPAMPath);
	size_t n = topoChecker.nTax();

#pragma omp parallel for
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t k = j + 1; k < n; ++k) {
				for (size_t l = k + 1; l < n; ++l) {
					std::array < size_t, 3 > counts = topoChecker.findMultiSPAMTopologyCounts(i, j, k, l);
					if (counts[0] + counts[1] + counts[2] > 0) {
						if (topoChecker.sameTopologyAsReference(i, j, k, l, counts)) {
#pragma omp atomic
							correct++;
						} else {
#pragma omp atomic
							wrong++;
						}
					}
				}
			}
		}
	}

	std::cout << "Finished the comparison! :-)\n";
	std::cout << "Correct topology multi-SpaM: " << correct << "\n";
	std::cout << "Wrong topology multi-SpaM: " << wrong << "\n";
	std::cout << "Quartet accuracy multi-SpaM: " << (double) correct / (correct + wrong) << "\n";
}


int main(int argc, char* argv[]) {
	if (argc != 4) {
		std::cout << "Usage: ./quartet_check path_to_fasta_file path_to_multispam_quartets path_to_reference_tree\n";
		return 1;
	}

	std::string fastaPath = argv[1];
	std::string quartetsPath = argv[2];
	std::string treePath = argv[3];
	evaluateQuartets(fastaPath, treePath, quartetsPath);
}
