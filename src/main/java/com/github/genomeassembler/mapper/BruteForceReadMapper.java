package com.github.genomeassembler.mapper;

import java.util.ArrayList;
import java.util.List;

public class BruteForceReadMapper implements ReadMapper {
    private final String genome;

    public BruteForceReadMapper(String genomeSequence) {
        this.genome = genomeSequence;
    }



    @Override
    public List<Integer> mapRead(String read, int mismatches) {
        List<Integer> startingPositions = new ArrayList<>();
        for (int i = 0; i < this.genome.length() - read.length() + 1; i++) {
            if (this.genome.substring(i, i + read.length()).equals(read)) {
                startingPositions.add(i);
            }
        }
        return startingPositions;
    }

    @Override
    public String getGenomeSequence() {
        return this.genome;
    }
}
