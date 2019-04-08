/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.arrays.illumina;
import org.apache.commons.io.IOUtils;
import picard.PicardException;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;


public class InfiniumEGTFile extends InfiniumDataFile {
    public static final String EXTENSION = "egt";

    private static final int VALID_GENTRAIN_DATA_TYPE = 9;
    private static final int INVALID_FILE_VERSION = 2;

    public int[] nAA;
    public int[] nAB;
    public int[] nBB;
    public float[] meanRAA;
    public float[] meanThetaAA;
    public float[] devRAA;
    public float[] devThetaAA;
    public float[] meanRAB;
    public float[] meanThetaAB;
    public float[] devRAB;
    public float[] devThetaAB;
    public float[] meanRBB;
    public float[] meanThetaBB;
    public float[] devRBB;
    public float[] devThetaBB;

    public float[] totalScore;

    public String[] rsNames;

    public String manifestName;

    public Map<String, Integer> rsNameToIndex;


    public InfiniumEGTFile(final File clusterFile) throws IOException {
        super(new DataInputStream(new FileInputStream(clusterFile)), false);
        parse();
    }

    private void parse() throws IOException {
        try {
            readHeaderData();
            readFileData();
        } finally {
            IOUtils.closeQuietly(stream);
        }
    }

    private void readFileData() throws IOException {
        final int version = parseInt();
        if (version > VALID_GENTRAIN_DATA_TYPE) {
            throw new IOException("Error. Cannot read file - unknown version " + version + " for gentrain data type");
        }

        manifestName = parseString();
        final int numCodes = parseInt();

        initializeArrays(numCodes);

        for (int i = 0; i < numCodes; i++) {
            nAA[i] = parseInt();
            nAB[i] = parseInt();
            nBB[i] = parseInt();
            devRAA[i] = parseFloat();
            devRAB[i] = parseFloat();
            devRBB[i] = parseFloat();
            meanRAA[i] = parseFloat();
            meanRAB[i] = parseFloat();
            meanRBB[i] = parseFloat();
            devThetaAA[i] = parseFloat();
            devThetaAB[i] = parseFloat();
            devThetaBB[i] = parseFloat();
            meanThetaAA[i] = parseFloat();
            meanThetaAB[i] = parseFloat();
            meanThetaBB[i] = parseFloat();

            // 15 unused floats
            skipFloats(15);
        }
        for (int i = 0; i < numCodes; i++) {
            // skip cSepScore
            skipFloat();

            totalScore[i] = parseFloat();

            // skip original score
            skipFloat();

            skipBoolean();
        }
        for (int i = 0; i < numCodes; i++) {
            skipString();
        }
        for (int i = 0; i < numCodes; i++) {
            rsNames[i] = parseString();
            if (rsNameToIndex.put(rsNames[i], i) != null) {
                throw new PicardException("Non-unique rsName '" + rsNames[i] + "' found in cluster file");
            }
        }
    }

    private void initializeArrays(int numCodes) {

        nAA = new int[numCodes];
        nAB = new int[numCodes];
        nBB = new int[numCodes];

        devRAA = new float[numCodes];
        devRAB = new float[numCodes];
        devRBB = new float[numCodes];
        meanRAA = new float[numCodes];
        meanRAB = new float[numCodes];
        meanRBB = new float[numCodes];
        devThetaAA = new float[numCodes];
        devThetaAB = new float[numCodes];
        devThetaBB = new float[numCodes];
        meanThetaAA = new float[numCodes];
        meanThetaAB = new float[numCodes];
        meanThetaBB = new float[numCodes];

        totalScore = new float[numCodes];
        rsNames = new String[numCodes];
        rsNameToIndex = new HashMap<>();
    }

    private void readHeaderData() throws IOException {
        setFileVersion(parseInt());
        // skip gcVersion
        skipString();
        // skip clusterVersion
        skipString();
        // skip callVersion
        skipString();
        // skip normalizationVersion
        skipString();
        // skip dataCreated
        skipString();
        // skip isWGT
        skipBoolean();

        if (getFileVersion() == INVALID_FILE_VERSION) {
            throw new IOException("Version '" + INVALID_FILE_VERSION + "' unsupported");
        }
        manifestName = parseString();
    }
}
