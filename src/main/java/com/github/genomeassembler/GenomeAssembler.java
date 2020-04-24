package com.github.genomeassembler;

import com.github.genomeassembler.mapper.ReadMapper;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GenomeAssembler {
    private final static String PROCESS_BEGIN_MSG = "Beginning processing of reference genome...";
    private final static String PROCESS_FINISH_MSG = "Finished processing of reference genome";
    private final static String PROCESSING_READS_MSG = "Retrieving reads from file...";

    private final ReadMapper referenceGenomeReadMapper;
    private final List<String> unmappedReads;
    private final Map<String, List<Integer>> exactMappedReads;


    public GenomeAssembler(String referenceGenomeSequence, List<String> reads) {
        // Process genome
        System.out.println("Length of genome: " + referenceGenomeSequence.length());
        System.out.println(PROCESS_BEGIN_MSG);
        this.referenceGenomeReadMapper = new ReadMapper(referenceGenomeSequence);
        System.out.println(PROCESS_FINISH_MSG);

        // Initialize reads
        this.unmappedReads = reads;
        this.exactMappedReads = new HashMap<>();
    }


    /**
     * Attempts to map every read with no tolerant mismatches.
     * If a read maps exactly, the read is removed from the unmappedReads list
     * to the mappedReads map
     * @return The number of reads (including duplicates) that were mapped
     *         to the reference genome
     */
    public int mapExactReads() {
        int count = 0;
        for (int i = this.unmappedReads.size() - 1; i >= 0; i--) {
            String read = this.unmappedReads.get(i);
            List<Integer> startingPositions =
                    this.referenceGenomeReadMapper.mapRead(read, 0);
            if (!startingPositions.isEmpty()) {
                this.unmappedReads.remove(i);
                exactMappedReads.put(read, startingPositions);
                count++;
            }
        }
        return count;
    }


    /**
     * Attempts to map remaining unmapped reads' reverse complements to the
     * reference genome with no tolerant mismatches.
     * If a read's reverse complement is mapped, the read is removed from
     * the unmapped reads, and its reverse complement is added to mapped
     * @return The number of reverse complementary reads (including duplicates)
     *         that were mapped to the reference genome
     */
    public int mapExactReverseComplementaryReads() {

    }









    /**
     * Retrieves the assembler's unmapped reads
     * @return A list containing all unmapped reads
     */
    public List<String> getUnmappedReads() {
        return this.unmappedReads;
    }


    /**
     * Retrieves the assembler's mapped reads
     * @return A map containing each mapped read and their respective starting
     *      positions in the reference genome
     */
    public Map<String, List<Integer>> getExactMappedReads() {
        return this.exactMappedReads;
    }






}
