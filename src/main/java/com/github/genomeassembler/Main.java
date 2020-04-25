package com.github.genomeassembler;

import com.github.genomeassembler.mapper.ReadMapper;
import com.github.genomeassembler.parser.FastaParser;
import com.github.genomeassembler.parser.FastqParser;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class Main {
    private final static String READS_PATH_ONE = "src/main/resources/";
    private final static String READS_PATH_TWO = "src/main/resources/";
    private final static String GENOME_PATH = "src/main/resources/M30.gff";

    public static void main(String[] args) {
        try {
            String genome = FastaParser.parseGffFile(GENOME_PATH);
            List<String> readsOne = FastqParser.getReads(READS_PATH_ONE);
            List<String> readsTwo = FastqParser.getReads(READS_PATH_TWO);

            List<String> combinedReads = new ArrayList<>(readsOne.size() + readsTwo.size());
            combinedReads.addAll(readsOne);
            combinedReads.addAll(readsTwo);
            GenomeAssembler genomeAssembler = new GenomeAssembler(genome, combinedReads);




        } catch (Exception e) {
            e.printStackTrace();
        }


        /**
         * Get contigs that map perfectly
         * For reads that don't map well, form contigs out of reads that have 1 - 3 by mapping each level
         * Finally, use de bruijn to get resulting contigs out of bad reads
         *
         * Map all contigs
         *
         */

    }
}
