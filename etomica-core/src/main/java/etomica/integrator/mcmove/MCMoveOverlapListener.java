/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.box.Box;
import etomica.math.numerical.AkimaSpline;
import etomica.util.IListener;

import java.util.Arrays;

public class MCMoveOverlapListener implements IListener<MCMoveEvent> {

    protected final MCMoveInsertDeleteBiased mcMove;
    protected double[][] sumInsert, sumDelete;
    protected long[] numInsert, numDelete;
    protected int minNumAtoms;
    protected double[] lnAlpha;
    protected double temperature;
    protected double[] ratios;
    protected final double[] y;
    protected final AkimaSpline spline;
    protected final int numAtomsLattice;
    protected double daDef, daSpan;
    protected int numAlpha;

    public MCMoveOverlapListener(MCMoveInsertDeleteBiased mcMove, int numAlpha, double daDef, int numAtomsLattice, double daSpan) {
        this.mcMove = mcMove;
        sumInsert = new double[0][0];
        sumDelete = new double[0][0];
        numInsert = new long[0];
        numDelete = new long[0];
        y = new double[numAlpha];
        spline = new AkimaSpline();
        ratios = new double[0];
        minNumAtoms = Integer.MAX_VALUE;
        this.numAtomsLattice = numAtomsLattice;
        this.daDef = daDef;
        this.daSpan = daSpan;
        this.numAlpha = numAlpha;
        lnAlpha = new double[numAlpha];
    }
    
    public void reset() {
        sumInsert = new double[0][0];
        sumDelete = new double[0][0];
        numInsert = new long[0];
        numDelete = new long[0];
        ratios = new double[0];
        minNumAtoms = Integer.MAX_VALUE;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    /**
     * Returns overlap sampling results (ratio of partition functions between
     * sequential number of atoms).  A ratio less than 1 means the free
     * energy biases the system towards having a vacancy.
     */
    public double[] getRatios() {
        if (minNumAtoms == Integer.MAX_VALUE) return ratios;
        if (ratios.length < (sumDelete.length-1) - minNumAtoms) {
            ratios = new double[(sumDelete.length-1) - minNumAtoms];
        }
        double[] z = new double[]{0};
        for (int na=minNumAtoms; na<sumDelete.length-1; na++) {
            int i = na - minNumAtoms;
            double da = daDef + Math.log(((double)(ratios.length-i))/(minNumAtoms+i+1));
            for (int j=0; j<numAlpha; j++) {
                // needs to be negative so y is increasing function of j
                double jLnAlpha = (da + daSpan*(j-(numAlpha-1)*0.5)*2.0/(numAlpha-1.0));
                lnAlpha[j] = jLnAlpha;
                y[j] = -Math.log((sumInsert[na][j]/numInsert[na]) / (sumDelete[na+1][j]/numDelete[na+1])) + jLnAlpha;
            }
            if (y[0] * y[y.length-1] > 0) {
//                System.out.println("failure "+na+" "+(na+1));
                if (!Double.isInfinite(y[0])) {
                    if (Math.abs(y[0]) < Math.abs(y[y.length-1])) {
                        double lnAlpha0 = (da + daSpan*(-(numAlpha-1)*0.5)*2.0/(numAlpha-1.0));
                        ratios[i] = Math.exp(lnAlpha0);
                    }
                    else {
                        double lnAlphaN = (da + daSpan*((numAlpha-1)-(numAlpha-1)*0.5)*2.0/(numAlpha-1.0));
                        ratios[i] = Math.exp(lnAlphaN);
                    }
                }
                else {
                    ratios[i] = Double.NaN;
                }
                continue;
            }
            spline.setInputData(y, lnAlpha);
            ratios[i] = Math.exp(spline.doInterpolation(z)[0]);
//            System.out.println("ratio: "+Math.log(ratios[i]));
        }
//        System.out.println(mcMove.species.getIndex()+" "+Arrays.toString(ratios));
        return ratios;
    }
    
    public double[] getHistogram() {
        if (minNumAtoms == Integer.MAX_VALUE) return new double[]{Double.NaN};
        double[] h = new double[sumDelete.length - minNumAtoms];
        double tot = 0;
        for (int na=minNumAtoms; na<sumDelete.length; na++) {
            int i = na - minNumAtoms;
            h[i] = numInsert[na] + numDelete[na];
            if (na == mcMove.minN || na == mcMove.maxN) {
                // numDelete[minN] = numInsert[maxN] = 0
                // we don't bother trying to jump of the ends, so double the other
                h[i] *= 2;
            }
            tot += h[i];
        }
        for (int i=0; i<h.length; i++) {
            h[i] /= tot;
        }
        return h;
    }
    
    public long[] getNumInsert() {
        if (minNumAtoms == Integer.MAX_VALUE) return new long[0];
        long[] h = new long[sumDelete.length - minNumAtoms];
        for (int na=minNumAtoms; na<sumDelete.length; na++) {
            int i = na - minNumAtoms;
            h[i] = numInsert[na];
        }
        return h;
    }
    
    public long[] getNumDelete() {
        if (minNumAtoms == Integer.MAX_VALUE) return new long[0];
        long[] h = new long[sumDelete.length - minNumAtoms];
        for (int na=minNumAtoms; na<sumDelete.length; na++) {
            int i = na - minNumAtoms;
            h[i] = numDelete[na];
        }
        return h;
    }
    
    public int getMinNumAtoms() {
        return minNumAtoms;
    }

    public void actionPerformed(MCMoveEvent event) {
        if (event.getMCMove() != mcMove) return;
        Box box = mcMove.getBox();
        int numAtoms = box.getNMolecules(mcMove.getSpecies());
        if (event instanceof MCMoveTrialFailedEvent) {
            // trial failed, but we were still here.  we need to increment our sums here
            // for the histogram.
            if (sumInsert.length < numAtoms+1) {
                sumInsert = Arrays.copyOf(sumInsert, numAtoms+1);
                numInsert = Arrays.copyOf(numInsert, numAtoms+1);
                sumDelete = Arrays.copyOf(sumDelete, numAtoms+1);
                numDelete = Arrays.copyOf(numDelete, numAtoms+1);
            }
        }
        else if (event instanceof MCMoveTrialInitiatedEvent) {
            // x = V/N*Math.exp(-beta*deltaU)
            double x = mcMove.getChi(temperature) * Math.exp(-mcMove.getLnBiasDiff());
            if (mcMove.lastMoveInsert()) {
                numAtoms--;
            }
            if (minNumAtoms > numAtoms) minNumAtoms = numAtoms;
            if (sumInsert.length < numAtoms+1) {
                sumInsert = Arrays.copyOf(sumInsert, numAtoms+1);
                numInsert = Arrays.copyOf(numInsert, numAtoms+1);
                sumDelete = Arrays.copyOf(sumDelete, numAtoms+1);
                numDelete = Arrays.copyOf(numDelete, numAtoms+1);
            }
            if (sumInsert[numAtoms] == null) {
                sumInsert[numAtoms] = new double[numAlpha];
                sumDelete[numAtoms] = new double[numAlpha];
            }
            if (mcMove.lastMoveInsert()) {
                double da = daDef + Math.log(((double)(numAtomsLattice-numAtoms))/(numAtoms+1));
                for (int i=0; i<numAlpha; i++) {
                    double iLnAlpha = (da + daSpan*(i-(numAlpha-1)*0.5)*2.0/(numAlpha-1.0));
                    sumInsert[numAtoms][i] += 1.0/(1 + Math.exp(iLnAlpha)/x);
                }
                numInsert[numAtoms]++;
            }
            else {
                double da = daDef + Math.log(((double)(numAtomsLattice-numAtoms+1))/numAtoms);
                for (int i=0; i<numAlpha; i++) {
                    double iLnAlpha = (da + daSpan*(i-(numAlpha-1)*0.5)*2.0/(numAlpha-1.0));
                    sumDelete[numAtoms][i] += 1.0/(Math.exp(iLnAlpha) + 1.0/x);
                }
                numDelete[numAtoms]++;
            }
        }
    }
}
