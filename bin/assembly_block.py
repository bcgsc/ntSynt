#!/usr/bin/env python3

'''
Represents an assembly block for an assembly.
Multiple assembly blocks comprise a synteny block
'''

class AssemblyBlock:
    "An assembly block for a given assembly. The AssemblyBlock objects per assembly make up a SyntenyBlock"
    def __init__(self, k):
        "Instantiate the AssemblyBlock"
        self.contig_id = None
        self.minimizers = []
        self.ori = None
        self.k = k

    def get_block_start(self):
        "Get the starting coordinate of the assembly block"
        return min(self.minimizers[0].position, self.minimizers[-1].position)

    def get_block_end(self):
        "Get the end coordinate of the assembly block"
        return max(self.minimizers[0].position, self.minimizers[-1].position) + self.k

    def get_block_length(self):
        "Get the length of the assembly block"
        return self.get_block_end() - self.get_block_start()

    def get_block_terminal_mx(self):
        "Return the terminal minimizer hashes for the assembly block"
        return self.contig_id, self.minimizers[0], self.minimizers[-1]

    def get_block_contig_start_end(self):
        "Return the main information about the block: contig, start, end"
        return self.contig_id, self.get_block_start(), self.get_block_end()

    def get_block_internal_mx_hashes(self):
        "Return the internal minimizer hashes for the assembly block"
        return [mx_pos.mx for mx_pos in self.minimizers[1:-1]]

    def __eq__(self, other):
        return (self.contig_id == other.contig_id) and (self.get_block_start() == other.get_block_start()) \
                and (self.get_block_end() == other.get_block_end())
