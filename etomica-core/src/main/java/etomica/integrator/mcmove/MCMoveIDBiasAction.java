/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.IAction;
import etomica.integrator.IntegratorMC;

/**
 * IAction which takes data from an MCMoveOverlapListener about acceptance
 * probabilities for insert/delete moves and uses them to update biasing of
 * the move.
 * 
 * @author Andrew Schultz
 */
public class MCMoveIDBiasAction implements IAction {
    protected final int maxVacancy;
    protected double bmu;
    protected final MCMoveOverlapListener mcMoveOverlapMeter;
    protected final MCMoveInsertDeleteBiased mcMoveID;
    protected final int numAtoms;
    protected final IntegratorMC integratorMC;
    protected int maxNumAtoms;
    protected double lastDefaultdADef;
    protected double pullFactor;
    protected final double ratioMaxN = 1, maxBiasMinN = 2;
    protected boolean playDumb = false;
    protected double fixedDaDef = Double.NaN;

    public MCMoveIDBiasAction(IntegratorMC integratorMC, MCMoveInsertDeleteBiased mcMoveID, int maxVacancy, double bmu,
                              MCMoveOverlapListener mcMoveOverlapMeter, int numAtoms) {
        this.integratorMC = integratorMC;
        this.mcMoveID = mcMoveID;
        this.maxVacancy = maxVacancy;
        this.bmu = bmu;
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.numAtoms = numAtoms;
        pullFactor = 1;
    }
    
    public void setFixedDaDef(double daDef) {
        fixedDaDef = daDef;
    }

    /**
     * With "play dumb" on, the bias will simply try to produce a flat
     * histogram based on previously-measured free energy differences
     * (no fitting or defect free energy will be used).
     */
    public void setPlayDumb(boolean playDumb) {
        this.playDumb = playDumb;
    }

    /**
     * Sets a nominal bias for states that have not yet been visited.  The bias given is
     * exp(-mu), so this should actually get something like mu/kT.
     */
    public void setMu(double bmu) {
        this.bmu = bmu;
        actionPerformed();
    }

    public void setNMaxReweight(int maxNumAtoms) {
        this.maxNumAtoms = maxNumAtoms;
    }

    public void setPullFactor(double pullFactor) {
        this.pullFactor = pullFactor;
    }

    public void setDefaultDaDef(double daDef) {
        lastDefaultdADef = daDef;
    }

    public void actionPerformed() {
        double[] ratios = mcMoveOverlapMeter.getRatios();
        double[] hist = mcMoveOverlapMeter.getHistogram();
        long[] numInsert = mcMoveOverlapMeter.getNumInsert();
        long[] numDelete = mcMoveOverlapMeter.getNumDelete();
        int n0 = mcMoveOverlapMeter.getMinNumAtoms();
        double daDef = fixedDaDef;
        if (Double.isNaN(daDef)) {
            double daDefTot = 0;
            double wTot = 0;
            for (int i=0; i<ratios.length; i++) {
                if (!Double.isNaN(ratios[i])) {
                    double iDaDef = Math.log(ratios[i]*((n0+i+1)/(ratios.length-i)));
                    double w = 1.0/(1.0/numInsert[i]+1.0/numDelete[i+1]);
                    daDefTot += iDaDef*w;
                    wTot += w;
                }
            }
            daDef = daDefTot/wTot;
        }
        double lnr = 0;
        double[] lnbias = new double[ratios.length+1];
        for (int i=0; i<ratios.length; i++) {
            lnbias[i] = lnr;
            if (maxNumAtoms>0 && n0+i+1>maxNumAtoms) {
                lnr += bmu;
                continue;
            }
            double lnratio = daDef + Math.log(((double)(numAtoms - (n0+i)))/(n0+i+1));
            if (playDumb) {
                lnratio = Math.log(ratios[i]);
            }
            if (!Double.isNaN(lnratio)) {
                if (-lnratio > bmu) {
                    // bias needed for flat histogram is greater than chemical potential
                    // only apply enough bias so that hist[N]/hist[i] ~= 0.25
                    double d = -lnratio - Math.log(ratioMaxN)/(1L<<(numAtoms - (n0+i)));
                    if (d < bmu) {
                        d = bmu;
                    }
                    lnr += d;
                }
                else {
                    // apply bias to make the histogram flat
                    lnr -= lnratio;
                    double offset = (lnbias[i] + bmu) - lnr;
                    if (offset > Math.log(maxBiasMinN)) {
                        // we're applying a lot of bias at low N.  apply only enough so
                        // that p(i+1)/p(i) = maxBiasMinN
                        offset = Math.log(maxBiasMinN);
                    }
                    lnr += offset;
                }
            }
            else {
                lnr += bmu;
            }
            if (!playDumb && hist[i]*hist[i+1] == 0 && hist[i]+hist[i+1] > 10000) {
                long ni = numInsert[i]+numDelete[i]+1;
                long nip1 = numInsert[i+1]+numDelete[i+1]+1;
                lnr += pullFactor*Math.log(((double)ni)/((double)nip1));
            }
//            else if (hist[i+1]<hist[i]) {
//                if (-Math.log(hist[i]/hist[i+1]) > Math.log(4)) {
//                    lnr += 0.25*Math.log(hist[i]/hist[i+1]);
//                }
//            }
//            else if (Math.abs(Math.log(hist[i]/hist[i+1])) > Math.log(10)) {
//                lnr += 0.1*Math.log(hist[i]/hist[i+1]);
//            }
        }
        lnbias[ratios.length] = lnr;
        int minN = numAtoms - maxVacancy;
        int maxN = numAtoms;
        double nominalLnBias = lnbias[lnbias.length-1];
        mcMoveID.setLnBias(numAtoms, 0);
        for (int na=numAtoms-1; na>=n0; na--) {
            mcMoveID.setLnBias(na, lnbias[na-n0]-nominalLnBias);
        }
        lnr = lnbias[0];
        if (Double.isNaN(daDef)) {
//            lastDefaultdADef -= 1;
            daDef = lastDefaultdADef;
        }
        int start = n0<numAtoms ? (n0-1) : (numAtoms-1);
        for (int na=start; na>=minN; na--) {
            double lnratio = daDef + Math.log(((double)(numAtoms-na))/(na+1));
            if (playDumb) {
                lnratio = bmu;
            }
            if (lnratio < bmu) {
                double offset = -(lnratio + bmu);
                if (offset < -Math.log(2)) offset = -Math.log(2);
                lnratio += offset;
                lnr += lnratio;
            }
            else {
                lnr += lnratio;
            }
            mcMoveID.setLnBias(na, lnr);
        }
        if (n0>numAtoms) return;

        for (int na = n0+lnbias.length; na<=maxN; na++) {
            double lastLnB = lnbias[lnbias.length-1]+(na-(n0+lnbias.length-1))*bmu;
            mcMoveID.setLnBias(na, lastLnB);
        }
    }
}
