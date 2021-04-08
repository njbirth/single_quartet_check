#pragma once

#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <sdsl/suffix_arrays.hpp>

#include "../genesis/lib/genesis/genesis.hpp"
#include "tree_information.hpp"

class TopologyChecker {
public:
	TopologyChecker(const std::string& fastaPath, const std::string& speciesTreePath, const std::string& multiSPAMPath);
	bool sameTopologyAsReference(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx, const std::array<size_t, 3>& counts);
	std::array<size_t, 3> findMultiSPAMTopologyCounts(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	size_t nTax();
private:
	size_t findReferenceTopology(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx);
	size_t rootIdx;
	std::vector<size_t> fastaIdxToRefIdx;
	TreeInformation informationReferenceTree;
	std::unordered_map<std::string, size_t> taxonToReferenceID;

	sdsl::csa_wt<> fmIndexMultiSPAM;
};
