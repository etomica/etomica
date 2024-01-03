/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.box.Box;
import etomica.data.meter.MeterStructureFactor;
import etomica.space.Vector;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

public class StructureFactorComponentWriter implements DataSinkBlockAveragerSFac.Sink, StructorFactorComponentExtractor.StructureFactorComponentSink, IDataSinkBlockAvg, Statefull {

    protected final Box box;
    protected final int minInterval;
    protected final int maxSteps;
    // we need to keep this separately (instead of using savedSteps[i].length) because we allocate savedAta[i][maxSteps]
    protected int[] numSavedSteps;
    protected float[][][][] savedData;
    protected int[] numSavedAvg;
    protected float[][] savedAvg;
    protected int[][] waveVectors;
    protected float constantAvg = Float.NaN;
    protected int grid = 10;

    public StructureFactorComponentWriter(MeterStructureFactor meterSFac, Box box, int minInterval, int maxSteps) {
        this.box = box;
        this.minInterval = minInterval;
        this.maxSteps = maxSteps;
        numSavedSteps = new int[0];
        savedData = new float[0][0][0][0];
        Vector L = box.getBoundary().getBoxSize();
        Vector[] wv = meterSFac.getWaveVectors();
        Vector wvL = box.getSpace().makeVector();
        waveVectors = new int[wv.length][wv[0].getD()];
        for (int j=0; j<wv.length; j++) {
            wvL.E(wv[j]);
            wvL.TE(L);
            wvL.TE(1/(2*Math.PI));
            for (int i=0; i<wvL.getD(); i++) {
                waveVectors[j][i] = (int)Math.round(wvL.getX(i));
            }
        }
        savedAvg = new float[0][];
        numSavedAvg = new int[0];
    }

    public void setConstantAverage(double avg) {
        this.constantAvg = (float)avg;
    }

    @Override
    public void putData(int interval, double[][] xy) {
        if (interval < minInterval) return;
        if (interval>=numSavedSteps.length) {
            numSavedSteps = Arrays.copyOf(numSavedSteps, interval+1);
            savedData = Arrays.copyOf(savedData, interval+1);
            savedData[interval] = new float[maxSteps][][];
        }
        if (maxSteps>0 && numSavedSteps[interval] >= maxSteps) return;
        savedData[interval][numSavedSteps[interval]] = new float[waveVectors.length][];
        int N = box.getLeafList().size();
        for (int i=0; i<xy.length; i++) {
            float beta = (float)(2*Math.sqrt((xy[i][0]*xy[i][0] + xy[i][1]*xy[i][1])*N));
            float theta = (float)Math.atan2(xy[i][0], xy[i][1]);
            savedData[interval][numSavedSteps[interval]][i] = new float[]{beta,theta};
        }
        numSavedSteps[interval]++;
    }

    public void writeFile(String filename) {
        try {
            float V = (float)box.getBoundary().volume();
            FileWriter fw = new FileWriter(filename);
            fw.write("{\"WV\": [");
            for (int j=0; j<waveVectors.length; j++) {
                if (j>0) fw.write(",");
                fw.write(" [");
                boolean firstD = true;
                for (int i=0; i<waveVectors[j].length; i++) {
                    if (!firstD) fw.write(",");
                    fw.write(""+waveVectors[j][i]);
                    firstD = false;
                }
                fw.write("]");
            }
            fw.write(" ],\n");
            if (savedAvg.length>0) {
                fw.write("\"avg\": [\n");
                for (int k=0; k< numSavedAvg.length; k++) {
                    if (numSavedAvg[k]*numSavedSteps[k] > 0 && numSavedAvg[k] != numSavedSteps[k]) throw new RuntimeException("Ooops, "+numSavedAvg[k]+" "+numSavedSteps[k]);
                    if (k>0) fw.write(",");
                    if (numSavedAvg[k]*numSavedSteps[k] == 0) {
                        fw.write("null");
                        continue;
                    }
                    fw.write("[");
                    for (int j=0; j<numSavedAvg[k]; j++) {
                        if (j>0) fw.write(",");
                        fw.write(String.format("%8.3e", savedAvg[k][j]));
                    }
                    fw.write("]\n");
                }
                fw.write("],\n");
            }
            else if (!Double.isNaN(constantAvg)) {
                fw.write("\"avg\": "+constantAvg+",\n");
            }

            fw.write("\"minmax\": [");
            for (int k=0; k<minInterval; k++) {
                if (k>0) fw.write(",");
                fw.write("null");
            }
            for (int k=minInterval; k<numSavedSteps.length; k++) {
                float vmin = Float.MAX_VALUE, vmax = -Float.MAX_VALUE;
                for (int j = 0; j < numSavedSteps[k]; j++) {
                    for (int ix = 0; ix<grid; ix++) {
                        float xx = ix/(float)grid;
                        for (int iy = 0; iy<grid; iy++) {
                            float yy = iy / (float) grid;
                            for (int iz = 0; iz<grid; iz++) {
                                float zz = iz / (float) grid;
                                float val = V*(savedAvg.length>0 ? savedAvg[k][j] : (Double.isNaN(constantAvg) ? 0 : constantAvg));
                                for (int i=0; i<waveVectors.length; i++) {
                                    int[] iwv = waveVectors[i];
                                    float beta = savedData[k][j][i][0];
                                    float theta = savedData[k][j][i][1];
                                    val += beta * Math.sin(2*Math.PI*(iwv[0]*xx + iwv[1]*yy + iwv[2]*zz) + theta);
                                }
                                if (val > vmax) vmax = val;
                                else if (val < vmin) vmin = val;
                            }
                        }
                    }
                }
                if (k>0) fw.write(",");
                fw.write(String.format("[ %8.3e, %8.3e ]", vmin/V, vmax/V));
            }
            fw.write("],\n");

            fw.write("\"sfac\": [");
            for (int k=0; k<minInterval; k++) {
                if (k>0) fw.write(",");
                fw.write("null");
            }
            for (int k=minInterval; k<numSavedSteps.length; k++) {
                fw.write(",[");
                for (int j=0; j<numSavedSteps[k]; j++) {
                    if (j>0) fw.write(",");
                    fw.write("[");
                    for (int i=0; i<waveVectors.length; i++) {
                        if (i>0) fw.write(",");
                        fw.write(String.format("[ %8.3e, %5.3f ]", savedData[k][j][i][0], savedData[k][j][i][1]));
                    }
                    fw.write("]");
                }
                fw.write("]\n");
            }

            fw.write("]}\n");
            fw.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public void putBlock(int interval, long step, double avg) {
        if (savedAvg.length <= interval) {
            savedAvg = Arrays.copyOf(savedAvg, interval+1);
            savedAvg[interval] = new float[maxSteps];
            numSavedAvg = Arrays.copyOf(numSavedAvg, interval+1);
        }
        if (numSavedAvg[interval] >= maxSteps) return;
        savedAvg[interval][numSavedAvg[interval]] = (float)avg;
        numSavedAvg[interval]++;
    }

    @Override
    public void putData(int idx, int interval, double[][] xyData) {
        // so ugly.  we take idx here because DataSourceCorrelation wants it
        putData(interval, xyData);
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        fw.write((numSavedSteps.length-1)+"\n");
        for (int i=0; i<numSavedSteps.length; i++) {
            if  (i>0) fw.write(" ");
            fw.write(numSavedSteps[i]+"");
        }
        fw.write("\n");
        for (int i=0; i<numSavedSteps.length; i++) {
            if (numSavedSteps[i]==0) continue;
            fw.write(i+" "+numSavedSteps[i]+"\n");
            for (int j=0; j<numSavedSteps[i]; j++) {
                for (int k=0; k<waveVectors.length; k++) {
                    if (k>0) fw.write(" ");
                    fw.write(savedData[i][j][k][0]+" "+savedData[i][j][k][1]);
                }
                fw.write("\n");
            }
        }
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        int maxInterval = Integer.parseInt(br.readLine());
        numSavedSteps = new int[maxInterval+1];
        String s = br.readLine();
        String[] bits = s.split(" ");
        for (int i=0; i<=maxInterval; i++) {
            numSavedSteps[i] = Integer.parseInt(bits[i]);
        }
        savedData = new float[maxInterval+1][][][];
        for (int i=0; i<numSavedSteps.length; i++) {
            if (numSavedSteps[i]==0) continue;
            bits = br.readLine().split(" ");
            int ii = Integer.parseInt(bits[0]);
            if (ii!=i) throw new RuntimeException("interval from file doesn't match");
            savedData[i] = new float[maxSteps][][];
            for (int j=0; j<numSavedSteps[i]; j++) {
                savedData[i][j] = new float[waveVectors.length][2];
                bits = br.readLine().split(" ");
                for (int k=0, l=0; k<waveVectors.length; k++, l+=2) {
                    savedData[i][j][k][0] = Float.parseFloat(bits[l]);
                    savedData[i][j][k][1] = Float.parseFloat(bits[l+1]);
                }
            }
        }
    }
}
