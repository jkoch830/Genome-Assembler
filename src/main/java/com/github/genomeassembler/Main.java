package com.github.genomeassembler;

import com.github.genomeassembler.parser.FastaParser;
import com.github.genomeassembler.parser.FastqParser;

import java.util.ArrayList;
import java.util.List;

public class Main {
    private final static String READS_PATH_ONE = "../832a.fastq";
    private final static String READS_PATH_TWO = "../832b.fastq";
    private final static String GENOME_PATH = "src/main/resources/M2282.gff";

    public static void main(String[] args) {
        try {
            String genome = FastaParser.parseGffFile(GENOME_PATH);
            List<String> readsOne = FastqParser.getReads(READS_PATH_ONE);
            List<String> readsTwo = FastqParser.getReads(READS_PATH_TWO);
            List<String> combinedReads = new ArrayList<>(readsOne.size() + readsTwo.size());
            combinedReads.addAll(readsOne);
            combinedReads.addAll(readsTwo);
            GenomeAssembler genomeAssembler = new GenomeAssembler(genome, combinedReads);
            AssemblerParameters parameters = new AssemblerParameters.Builder().
                    mismatchToleranceLowerBound(0).
                    mismatchToleranceHigherBound(7).
                    mismatchToleranceStep(1).
                    minContigOutputLength(250).
                    kmerLength(30).
                    build();
            genomeAssembler.setAssemblerParameters(parameters);
            genomeAssembler.assemble();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
