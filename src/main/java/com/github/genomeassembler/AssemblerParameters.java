package com.github.genomeassembler;

import java.util.Objects;

/**
 * Holds the genome assembler's settings
 */
public class AssemblerParameters {
    private final int requiredContigOverlap;
    private final int minContigOutputLength;
    private final int mismatchToleranceLowerBound;
    private final int mismatchToleranceHigherBound;
    private final int mismatchToleranceStep;
    private final int kmerLength;

    /**
     * Constructor following builder method
     */
    public static class Builder {

        // Optional parameters - initialized to default values
        private int requiredContigOverlap = 0;
        private int minContigOutputLength = 1000;
        private int mismatchToleranceLowerBound = 0;
        private int mismatchToleranceHigherBound = 15;
        private int mismatchToleranceStep = 3;
        private int kmerLength = 30;

        public Builder requiredContigOverlap(int val) {
            requiredContigOverlap = val;
            return this;
        }

        public Builder minContigOutputLength(int val) {
            minContigOutputLength = val;
            return this;
        }

        public Builder mismatchToleranceLowerBound(int val) {
            mismatchToleranceLowerBound = val;
            return this;
        }

        public Builder mismatchToleranceHigherBound(int val) {
            mismatchToleranceHigherBound = val;
            return this;
        }

        public Builder mismatchToleranceStep(int val) {
            mismatchToleranceStep = val;
            return this;
        }

        public Builder kmerLength(int val) {
            kmerLength = val;
            return this;
        }

        public AssemblerParameters build() {
            return new AssemblerParameters(this);
        }
    }

    private AssemblerParameters(Builder builder) {
        requiredContigOverlap = builder.requiredContigOverlap;
        minContigOutputLength = builder.minContigOutputLength;
        mismatchToleranceLowerBound = builder.mismatchToleranceLowerBound;
        mismatchToleranceHigherBound = builder.mismatchToleranceHigherBound;
        mismatchToleranceStep = builder.mismatchToleranceStep;
        kmerLength = builder.kmerLength;
    }

    public int getRequiredContigOverlap() {
        return requiredContigOverlap;
    }

    public int getMinContigOutputLength() {
        return minContigOutputLength;
    }

    public int getMismatchToleranceLowerBound() {
        return mismatchToleranceLowerBound;
    }

    public int getMismatchToleranceHigherBound() {
        return mismatchToleranceHigherBound;
    }

    public int getMismatchToleranceStep() {
        return mismatchToleranceStep;
    }

    public int getKmerLength() {
        return kmerLength;
    }


}
