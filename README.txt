
Compilation:
mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make

Usage: ./quartet_check path_to_fasta_file path_to_multispam_quartets path_to_reference_tree

~~Unfortunately, the code still has a dependency on the sdsl-lite library (https://github.com/simongog/sdsl-lite) which needs to be installed manually.~~

EXAMPLE:
./quartet_check example_data/59_unaligned.fasta example_data/59_multiSPAM_quartets.txt example_data/59.phy.raxml.bestTree
