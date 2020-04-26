package com.github.genomeassembler.mapper;

import java.util.List;

/**
 * Interface of a read-mapper. A read-mapper pre-processes a genome to
 * allow efficient mapping of reads to the genome.
 */
public interface ReadMapper {

    /**
     * Maps a read to the genome
     * @param read The read being mapped
     * @param mismatches The number of tolerant mismatches
     * @return All starting positions in the genome that read occurs
     */
    List<Integer> mapRead(String read, int mismatches);

    /**
     * Gets the genome sequence
     * @return The genome sequence
     */
    String getGenomeSequence();
}
