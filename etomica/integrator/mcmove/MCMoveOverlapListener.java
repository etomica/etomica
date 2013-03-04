package etomica.integrator.mcmove;

import etomica.api.IBox;
import etomica.util.IEvent;
import etomica.util.IListener;
import etomica.util.numerical.AkimaSpline;

public class MCMoveOverlapListener implements IListener {

    protected final MCMoveInsertDeleteBiased mcMove;
    protected double[][] sumInsert, sumDelete;
    protected double[] numInsert, numDelete;
    protected int minNumAtoms;
    protected double[] alpha, lnAlpha;
    protected double temperature;
    protected double[] ratios;
    protected double[] y;
    protected final AkimaSpline spline;
    
    public MCMoveOverlapListener(MCMoveInsertDeleteBiased mcMove, int numAlpha, double alphaCenter, double alphaSpan) {
        this.mcMove = mcMove;
        sumInsert = new double[0][0];
        sumDelete = new double[0][0];
        numInsert = new double[0];
        numDelete = new double[0];
        alpha = new double[numAlpha];
        lnAlpha = new double[numAlpha];
        for (int i=0; i<numAlpha; i++) {
            double dx = numAlpha==1 ? 0 : (i - (numAlpha-1)*0.5)/(numAlpha-1)*2.0;
            alpha[i] = alphaCenter * Math.pow(alphaSpan, dx);
            lnAlpha[i] = Math.log(alpha[i]);
        }
        y = new double[numAlpha];
        spline = new AkimaSpline();
        ratios = new double[0];
        minNumAtoms = Integer.MAX_VALUE;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    public double[] getRatios() {
        if (ratios.length < (sumDelete.length-1) - minNumAtoms) {
            ratios = new double[(sumDelete.length-1) - minNumAtoms];
        }
        double[] z = new double[]{0};
        if (minNumAtoms == Integer.MAX_VALUE) return ratios;
        for (int na=minNumAtoms; na<sumDelete.length-1; na++) {
            int i = na - minNumAtoms;
            for (int j=0; j<alpha.length; j++) {
                double jAlpha = alpha[j];
                // needs to be negative so y is increasing function of j
                y[j] = -Math.log((sumInsert[na][j]/numInsert[na]) / (sumDelete[na+1][j]/numDelete[na+1]) / jAlpha);
            }
            if (y[0] * y[y.length-1] > 0) {
//                System.out.println("failure "+na+" "+(na+1));
                if (!Double.isInfinite(y[0])) {
                    if (Math.abs(y[0]) < Math.abs(y[y.length-1])) {
                        ratios[i] = alpha[0];
                    }
                    else {
                        int last = alpha.length-1;
                        ratios[i] = alpha[last-1];
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
        return ratios;
    }
    
    public double[] getHistogram() {
        double[] h = new double[sumDelete.length - minNumAtoms];
        double tot = 0;
        for (int na=minNumAtoms; na<sumDelete.length; na++) {
            int i = na - minNumAtoms;
            h[i] = numInsert[na] + numDelete[na];
            tot += h[i];
        }
        for (int i=0; i<h.length; i++) {
            h[i] /= tot;
        }
        return h;
    }
    
    public int getMinNumAtoms() {
        return minNumAtoms;
    }

    public void actionPerformed(IEvent event) {
        if (event instanceof MCMoveTrialFailedEvent) {
            if (((MCMoveEvent)event).getMCMove() != mcMove) return;
            double w = 1.0;
            IBox box = mcMove.getBox();
            int numAtoms = box.getLeafList().getAtomCount();
            if (sumInsert.length < numAtoms+1) {
                sumInsert = (double[][])etomica.util.Arrays.resizeArray(sumInsert, numAtoms+1);
                numInsert = etomica.util.Arrays.resizeArray(numInsert, numAtoms+1);
                sumDelete = (double[][])etomica.util.Arrays.resizeArray(sumDelete, numAtoms+1);
                numDelete = etomica.util.Arrays.resizeArray(numDelete, numAtoms+1);
            }
            if (mcMove.lastMoveInsert()) {
                numInsert[numAtoms] += w;
            }
            else {
                numDelete[numAtoms] += w;
            }
        }
        else if (event instanceof MCMoveTrialInitiatedEvent) {
            if (((MCMoveEvent)event).getMCMove() != mcMove) return;
            double w = 1.0;
            IBox box = mcMove.getBox();
            int numAtoms = box.getLeafList().getAtomCount();
            double x = mcMove.getA();
            x *= Math.exp(-mcMove.getLnBiasDiff() + mcMove.getB()/temperature);
            if (mcMove.lastMoveInsert()) {
                numAtoms--;
            }
            if (minNumAtoms > numAtoms) minNumAtoms = numAtoms;
            if (sumInsert.length < numAtoms+1) {
                sumInsert = (double[][])etomica.util.Arrays.resizeArray(sumInsert, numAtoms+1);
                numInsert = etomica.util.Arrays.resizeArray(numInsert, numAtoms+1);
                sumDelete = (double[][])etomica.util.Arrays.resizeArray(sumDelete, numAtoms+1);
                numDelete = etomica.util.Arrays.resizeArray(numDelete, numAtoms+1);
            }
            if (sumInsert[numAtoms] == null) {
                sumInsert[numAtoms] = new double[alpha.length];
                sumDelete[numAtoms] = new double[alpha.length];
            }
            if (mcMove.lastMoveInsert()) {
                for (int i=0; i<alpha.length; i++) {
                    double iAlpha = alpha[i];
                    sumInsert[numAtoms][i] += w/(1 + iAlpha/x);
                }
                numInsert[numAtoms] += w;
            }
            else {
                for (int i=0; i<alpha.length; i++) {
                    double iAlpha = alpha[i];
                    sumDelete[numAtoms][i] += w/(iAlpha + 1.0/x);
                }
                numDelete[numAtoms] += w;
            }
        }
    }
}
