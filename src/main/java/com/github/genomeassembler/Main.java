package com.github.genomeassembler;

import com.github.genomeassembler.mapper.ReadMapper;
import com.github.genomeassembler.parser.FastaParser;
import com.github.genomeassembler.parser.FastqParser;

import java.util.List;
import java.util.Map;

public class Main {
    private static String READS_PATH;
    private static String GENOME_PATH;
    public static void main(String args[]) {
        String genome = FastaParser.parseFasta(GENOME_PATH);
        List<String> reads = FastqParser.getReads(READS_PATH);
        ReadMapper readMapper = new ReadMapper(genome);

        // Finds exact matches of reads against genome
        Map<String, List<Integer>> mappedReads = readMapper.mapReads(reads, 0);



        /**
         * First align exact reads (form contigs)
         *      Try each reverse complement if it doesn't match
         * Then try to align reads with up to d mismatches
         * De novo assembly of rest of reads
         * Align contigs
         *
         *
         *
         * Every Read will be an int[], int[0] is length
         *
         */

    }
}
