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

size_t TopologyChecker::nTax() {
	return fastaIdxToRefIdx.size();
}

TopologyChecker::TopologyChecker(const std::string& fastaPath, const std::string& speciesTreePath,
		const std::string& multiSPAMPath) {
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

	std::vector<std::string> fastaLabels = readFastaLabels(fastaPath);

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

	if (!multiSPAMPath.empty()) {
		sdsl::construct(fmIndexMultiSPAM, multiSPAMPath, 1);
	}
}

size_t TopologyChecker::findReferenceTopology(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
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

	if (informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			> informationReferenceTree.distanceInEdges(lca_uw, lca_vz)
			&& informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
					> informationReferenceTree.distanceInEdges(lca_uz, lca_vw)) {
		return 0; // ab|cd
	} else if (informationReferenceTree.distanceInEdges(lca_uw, lca_vz)
			> informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			&& informationReferenceTree.distanceInEdges(lca_uw, lca_vz)
					> informationReferenceTree.distanceInEdges(lca_uz, lca_vw)) {
		return 1; // ac|bd
	} else if (informationReferenceTree.distanceInEdges(lca_uz, lca_vw)
			> informationReferenceTree.distanceInEdges(lca_uv, lca_wz)
			&& informationReferenceTree.distanceInEdges(lca_uz, lca_vw)
					> informationReferenceTree.distanceInEdges(lca_uw, lca_vz)) {
		return 2; // ad|bc
	} else {
		return 3; // star topology
	}
}

bool TopologyChecker::sameTopologyAsReference(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx,
		const std::array<size_t, 3>& counts) {
	size_t refTopology = findReferenceTopology(aIdx, bIdx, cIdx, dIdx);
	size_t q1, q2, q3;
	if (refTopology == 0) { // ab|cd
		q1 = 0;
		q2 = 1;
		q3 = 2;
	} else if (refTopology == 1) { // ac|bd
		q1 = 1;
		q2 = 0;
		q3 = 2;
	} else if (refTopology == 2) { // ad|bc
		q1 = 2;
		q2 = 1;
		q3 = 0;
	} else {
		throw std::runtime_error("The quartet has star-topology in the reference tree");
	}
	if (counts[q1] < counts[q2] || counts[q1] < counts[q3]) {
		return false;
	} else {
		return true;
	}
}

std::array<size_t, 3> TopologyChecker::findMultiSPAMTopologyCounts(size_t aIdx, size_t bIdx, size_t cIdx, size_t dIdx) {
	std::array < size_t, 3> counts = {0,0,0};
	std::string abcd = std::to_string(aIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(dIdx) + ":";
	std::string abdc = std::to_string(aIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(cIdx) + ":";
	std::string acbd = std::to_string(aIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(dIdx) + ":";
	std::string acdb = std::to_string(aIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(bIdx) + ":";
	std::string adbc = std::to_string(aIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(cIdx) + ":";
	std::string adcb = std::to_string(aIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(bIdx) + ":";

	std::string bacd = std::to_string(bIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(dIdx) + ":";
	std::string badc = std::to_string(bIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(cIdx) + ":";
	std::string bcad = std::to_string(bIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(dIdx) + ":";
	std::string bcda = std::to_string(bIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(aIdx) + ":";
	std::string bdac = std::to_string(bIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(cIdx) + ":";
	std::string bdca = std::to_string(bIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(aIdx) + ":";

	std::string cabd = std::to_string(cIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(dIdx) + ":";
	std::string cadb = std::to_string(cIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(bIdx) + ":";
	std::string cbad = std::to_string(cIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(dIdx) + ":";
	std::string cbda = std::to_string(cIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(dIdx) + "," + std::to_string(aIdx) + ":";
	std::string cdab = std::to_string(cIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(bIdx) + ":";
	std::string cdba = std::to_string(cIdx) + "," + std::to_string(dIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(aIdx) + ":";

	std::string dabc = std::to_string(dIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(cIdx) + ":";
	std::string dacb = std::to_string(dIdx) + "," + std::to_string(aIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(bIdx) + ":";
	std::string dbac = std::to_string(dIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(cIdx) + ":";
	std::string dbca = std::to_string(dIdx) + "," + std::to_string(bIdx) + "|" + std::to_string(cIdx) + "," + std::to_string(aIdx) + ":";
	std::string dcab = std::to_string(dIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(aIdx) + "," + std::to_string(bIdx) + ":";
	std::string dcba = std::to_string(dIdx) + "," + std::to_string(cIdx) + "|" + std::to_string(bIdx) + "," + std::to_string(aIdx) + ":";

	// ab|cd topology:
	counts[0] += sdsl::count(fmIndexMultiSPAM, abcd.begin(), abcd.end());
	counts[0] += sdsl::count(fmIndexMultiSPAM, abdc.begin(), abdc.end());
	counts[0] += sdsl::count(fmIndexMultiSPAM, bacd.begin(), bacd.end());
	counts[0] += sdsl::count(fmIndexMultiSPAM, badc.begin(), badc.end());
	counts[0] += sdsl::count(fmIndexMultiSPAM, cdab.begin(), cdab.end());
	counts[0] += sdsl::count(fmIndexMultiSPAM, cdba.begin(), cdba.end());
	counts[0] += sdsl::count(fmIndexMultiSPAM, dcab.begin(), dcab.end());
	counts[0] += sdsl::count(fmIndexMultiSPAM, dcba.begin(), dcba.end());
	// ac|bd topology:
	counts[1] += sdsl::count(fmIndexMultiSPAM, acbd.begin(), acbd.end());
	counts[1] += sdsl::count(fmIndexMultiSPAM, acdb.begin(), acdb.end());
	counts[1] += sdsl::count(fmIndexMultiSPAM, cabd.begin(), cabd.end());
	counts[1] += sdsl::count(fmIndexMultiSPAM, cadb.begin(), cadb.end());
	counts[1] += sdsl::count(fmIndexMultiSPAM, bdac.begin(), bdac.end());
	counts[1] += sdsl::count(fmIndexMultiSPAM, bdca.begin(), bdca.end());
	counts[1] += sdsl::count(fmIndexMultiSPAM, dbac.begin(), dbac.end());
	counts[1] += sdsl::count(fmIndexMultiSPAM, dbca.begin(), dbca.end());
	// ad|bc topology:
	counts[2] += sdsl::count(fmIndexMultiSPAM, adbc.begin(), adbc.end());
	counts[2] += sdsl::count(fmIndexMultiSPAM, adcb.begin(), adcb.end());
	counts[2] += sdsl::count(fmIndexMultiSPAM, dabc.begin(), dabc.end());
	counts[2] += sdsl::count(fmIndexMultiSPAM, dacb.begin(), dacb.end());
	counts[2] += sdsl::count(fmIndexMultiSPAM, bcad.begin(), bcad.end());
	counts[2] += sdsl::count(fmIndexMultiSPAM, bcda.begin(), bcda.end());
	counts[2] += sdsl::count(fmIndexMultiSPAM, cbda.begin(), cbda.end());
	counts[2] += sdsl::count(fmIndexMultiSPAM, cbad.begin(), cbad.end());

	return counts;
}
