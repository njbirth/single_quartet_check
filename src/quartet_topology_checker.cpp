#include "quartet_topology_checker.hpp"

#include <fstream>

std::vector<std::string> readFastaLabels(const std::string& fastaPath) {
	std::vector<std::string> res;
	std::ifstream infile(fastaPath);
	while (infile.good()) {
		std::string line;
		infile >> line;
		if (line[0] == '>') {
			res.push_back(line.substr(1, line.size()));
		}
	}
	infile.close();
	return res;
}

int TopologyChecker::labelToInt(std::string label) {
	for(int i = 0; i < fastaLabels.size(); i++) {
		if(fastaLabels[i].compare(label) == 0)
			return i;
	}

	return -1;
}

TopologyChecker::TopologyChecker(const std::string& fastaPath, const std::string& speciesTreePath) {
	// read the species tree; including Sebastians quickfix
	genesis::tree::Tree speciesTree = genesis::tree::DefaultTreeNewickReader().from_file(speciesTreePath);
	if (speciesTree.root_node().rank() == 1) { // unroot the tree if it is not already unrooted, code from Master student
		std::string newick = genesis::tree::DefaultTreeNewickWriter().to_string(speciesTree);
		int c = 0;
		size_t a = 0;
		size_t b = newick.size();
		size_t d = newick.size();
		for (size_t i = 0; i < newick.size() and b == newick.size(); ++i) {
			if (newick[i] == '(')
				c++;
			if (newick[i] == ')')
				c--;
			if (c == 2 and a == 0)
				a = i;
			if (c == 1 and a > 0) {
				b = i;
				d = b;
				if (newick[i + 1] == ':') {
					d++;
					while (newick[d] != ',' and newick[d] != ')')
						d++;
				}
			}
		}
		std::string newick2 = newick.substr(0, a) + newick.substr(a + 1, b - a - 1)
				+ newick.substr(d, newick.size() - d);
		speciesTree = genesis::tree::DefaultTreeNewickReader().from_string(newick2);
	}
	rootIdx = speciesTree.root_node().index();

	fastaLabels = readFastaLabels(fastaPath);

	fastaIdxToRefIdx.resize(fastaLabels.size());
	for (size_t i = 0; i < speciesTree.node_count(); ++i) {
		if (speciesTree.node_at(i).is_leaf()) {
			std::string leafName = speciesTree.node_at(i).data<genesis::tree::DefaultNodeData>().name;
			taxonToReferenceID[leafName] = i;
			bool found = false;
			for (size_t j = 0; j < fastaLabels.size(); ++j) {
				if (fastaLabels[j] == leafName) {
					fastaIdxToRefIdx[j] = i;
					found = true;
					break;
				}
			}
			if (!found) {
				throw std::runtime_error("Could not find any taxon named " + leafName + " in the FASTA file");
			}
		}
	}

	informationReferenceTree.init(speciesTree);
}

// Returns true if a,b|c,d is the correct topology according to the reference tree
bool TopologyChecker::sameTopologyAsReference(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	size_t uIdx = fastaIdxToRefIdx[aIdx];
	size_t vIdx = fastaIdxToRefIdx[bIdx];
	size_t wIdx = fastaIdxToRefIdx[cIdx];
	size_t zIdx = fastaIdxToRefIdx[dIdx];

	size_t lca_uv = informationReferenceTree.lowestCommonAncestorIdx(uIdx, vIdx, rootIdx);
	size_t lca_uw = informationReferenceTree.lowestCommonAncestorIdx(uIdx, wIdx, rootIdx);
	size_t lca_uz = informationReferenceTree.lowestCommonAncestorIdx(uIdx, zIdx, rootIdx);
	size_t lca_vw = informationReferenceTree.lowestCommonAncestorIdx(vIdx, wIdx, rootIdx);
	size_t lca_vz = informationReferenceTree.lowestCommonAncestorIdx(vIdx, zIdx, rootIdx);
	size_t lca_wz = informationReferenceTree.lowestCommonAncestorIdx(wIdx, zIdx, rootIdx);

	return informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			> informationReferenceTree.distanceInEdges(lca_uw, lca_vz)
		&& informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			> informationReferenceTree.distanceInEdges(lca_uz, lca_vw);
}