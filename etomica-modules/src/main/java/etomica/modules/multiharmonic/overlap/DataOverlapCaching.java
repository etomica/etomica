/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.overlap;

import java.util.HashMap;

import etomica.overlap.AlphaSource;
import etomica.overlap.DataOverlap;
import etomica.overlap.IntegratorOverlap;

/**
 * DataOverlap subclass that caches return values.
 * 
 * @author Andrew Schultz
 */
public class DataOverlapCaching extends DataOverlap {

    protected long[] lastStepCount = new long[4];
    protected final double[] lastResult = new double[2];
    protected final IntegratorOverlap integrator;
    protected final HashMap<String,double[]>[] seenParams;
    protected final long[][] lastStepPop;
    protected double[][][][] seenArrays;
    
    public DataOverlapCaching(DataSourceOverlapAvg refAvg, DataSourceOverlapAvg targetAvg, AlphaSource alphaSource, IntegratorOverlap integrator) {
        super(refAvg, targetAvg, alphaSource);
        this.integrator = integrator;
        seenParams = new HashMap[3];
        seenParams[0] = new HashMap<String,double[]>();
        seenParams[1] = new HashMap<String,double[]>();
        seenParams[2] = new HashMap<String,double[]>();
        lastStepCount[0] = -1;
        lastStepCount[1] = -1;
        lastStepCount[2] = -1;
        lastStepCount[3] = -1;
        lastStepPop = new long[][]{{-1,-1},{-1,-1},{-1,-1}};
        seenArrays = new double[3][2][4][0];
    }

    public double[] getOverlapAverageAndError() {
        long newStepCount = integrator.getStepCount();
        if (lastStepCount[3] == newStepCount) return lastResult;
        double[] superResult = super.getOverlapAverageAndError();
        lastResult[0] = superResult[0];
        lastResult[1] = superResult[1];
        lastStepCount[3] = newStepCount;
        return lastResult;
    }

    protected double[] getAverageAndError(int which, double iAlpha, boolean doLog) {
        long newStepCount = -1;
        if (which == RATIO) {
            newStepCount = integrator.getIntegrators()[0].getStepCount() + integrator.getIntegrators()[1].getStepCount();
        }
        else {
            newStepCount = integrator.getIntegrators()[which].getStepCount();
        }
        if (lastStepCount[which] != newStepCount) {
            seenParams[which].clear();
        }
        String hashStr = which+"_"+iAlpha+"_"+doLog;
        double[] result = seenParams[which].get(hashStr);
        if (result != null) return result;
        double[] superResult = super.getAverageAndError(which, iAlpha, doLog);
        result = new double[]{superResult[0],superResult[1]};
        seenParams[which].put(hashStr, result);
        lastStepCount[which] = newStepCount;
        return result;
    }

    protected void populateArrays(int which, boolean doLog) {
        long newStepCount = -1;
        int logNum = doLog ? 0 : 1;
        if (which == RATIO) {
            newStepCount = integrator.getIntegrators()[0].getStepCount() + integrator.getIntegrators()[1].getStepCount();
        }
        else {
            newStepCount = integrator.getIntegrators()[which].getStepCount();
        }
        int numAlpha = alphaSource.getNumAlpha();
        if (seenArrays[0][0][0].length != numAlpha) {
            seenArrays = new double[3][2][4][numAlpha];
        }
        if (lastStepPop[which][logNum] != newStepCount) {
            super.populateArrays(which, doLog);
            System.arraycopy(lnAlpha, 0, seenArrays[which][logNum][0], 0, numAlpha);
            System.arraycopy(err, 0, seenArrays[which][logNum][1], 0, numAlpha);
            System.arraycopy(lnAlphaDiff, 0, seenArrays[which][logNum][2], 0, numAlpha);
            System.arraycopy(lnRatio, 0, seenArrays[which][logNum][3], 0, numAlpha);
            
            lastStepPop[which][doLog?0:1] = newStepCount;
        }
        else {
            System.arraycopy(seenArrays[which][logNum][0], 0, lnAlpha, 0, numAlpha);
            System.arraycopy(seenArrays[which][logNum][1], 0, err, 0, numAlpha);
            System.arraycopy(seenArrays[which][logNum][2], 0, lnAlphaDiff, 0, numAlpha);
            System.arraycopy(seenArrays[which][logNum][3], 0, lnRatio, 0, numAlpha);
        }
    }
}
