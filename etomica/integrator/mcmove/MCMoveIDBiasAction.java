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
    private final int maxDN;
    private final int fixedN;
    private double mu;
    private final MCMoveOverlapListener mcMoveOverlapMeter;
    protected final MCMoveInsertDeleteBiased mcMoveID;
    private final int numAtoms;
    private final double temperature;
    protected final IntegratorMC integratorMC;

    public MCMoveIDBiasAction(IntegratorMC integratorMC, MCMoveInsertDeleteBiased mcMoveID, int maxDN, int fixedN, double mu,
            MCMoveOverlapListener mcMoveOverlapMeter, int numAtoms, double temperature) {
        this.integratorMC = integratorMC;
        this.mcMoveID = mcMoveID;
        this.maxDN = maxDN;
        this.fixedN = fixedN;
        this.mu = mu;
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.numAtoms = numAtoms;
        this.temperature = temperature;
    }

    /**
     * Sets a nominal bias for states that have not yet been visited.  The bias given is
     * exp(-mu), so this should actually get something like mu/kT.
     */
    public void setMu(double mu) {
        this.mu = mu;
        actionPerformed();
    }

    public void actionPerformed() {
        double[] ratios = mcMoveOverlapMeter.getRatios();
        if (ratios.length == 0) return;
        double[] hist = mcMoveOverlapMeter.getHistogram();
        double lnr = 0;
        int n0 = mcMoveOverlapMeter.getMinNumAtoms();
        double[] lnbias = new double[ratios.length+1];
        double[] lnbias2 = new double[ratios.length+1];
        for (int i=0; i<ratios.length; i++) {
            lnbias[i] = -lnr;
            if (!Double.isNaN(ratios[i])) {
                lnr += Math.log(ratios[i]);
            }
            else {
                lnr -= mu/temperature;
            }
        }
        lnbias[ratios.length] = -lnr;
        int minN = (fixedN < numAtoms ? fixedN : numAtoms) - maxDN - 1;
        int maxN = (fixedN > numAtoms ? fixedN : numAtoms) + maxDN + 1;
        double nominalLnBias = lnbias[numAtoms-n0];
        mcMoveID.setLnBias(n0-1, (lnbias[0]-nominalLnBias) - mu/temperature + 10);
        for (int i=0; i<lnbias.length; i++) {
            lnbias[i] -= nominalLnBias;
            lnbias2[i] = lnbias[i];
            if (hist[i]*hist[hist.length-1] != 0) lnbias2[i] -= 10*Math.log(hist[i]/hist[numAtoms-n0]);
            int na = n0 + i;
            mcMoveID.setLnBias(na, lnbias2[i]);
        }
        for (int na=minN; na<n0; na++) {
            mcMoveID.setLnBias(na, lnbias2[0] - (mu/temperature)*(n0-na));
        }

        for (int na = n0+lnbias.length; na<=maxN; na++) {
            double lastLnB = lnbias2[lnbias2.length-1]+(na-(n0+lnbias.length-1))*(mu/temperature);
            mcMoveID.setLnBias(na, lastLnB);
        }
        System.out.print("\n");
    }
}