#pragma once

#include <vector>
#include <array>
#include <string>
#include <unordered_map>

#include "../genesis/lib/genesis/genesis.hpp"
#include "tree_information.hpp"

class TopologyChecker {
public:
	TopologyChecker(const std::string& fastaPath, const std::string& speciesTreePath);
	bool sameTopologyAsReference(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	int labelToInt(std::string label);
private:
	size_t rootIdx;
	std::vector<size_t> fastaIdxToRefIdx;
	TreeInformation informationReferenceTree;
	std::unordered_map<std::string, size_t> taxonToReferenceID;
	std::vector<std::string> fastaLabels;
};
