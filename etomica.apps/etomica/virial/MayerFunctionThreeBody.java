package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.virial.cluster.VirialDiagrams;

/**
 * This class acts as a non-additive Mayer function for molecules that
 * interact with only a 3-body potential.
 * 
 * This class handles bookkeeping and relies on a subclass to know how to
 * actually calculate the potential.
 *
 * @author Andrew Schultz
 */
public abstract class MayerFunctionThreeBody implements MayerFunctionNonAdditive {

    protected double[] lastValue;
    protected int totalMolecules;

    /*
     * We assume that this will be called for all triplets before it is called
     * for anything else.  We'll calculate the higher order energies by summing
     * the triplet contributions
     */
    public double f(IMoleculeList molecules, double[] r2, double beta) {
        int nMolecules = molecules.getMoleculeCount();
        if (nMolecules<3) throw new RuntimeException("need at least 3");
        if (nMolecules==3) {
            int id0 = molecules.getMolecule(0).getIndex();
            int id1 = molecules.getMolecule(1).getIndex();
            int id2 = molecules.getMolecule(2).getIndex();
            int tripletID = VirialDiagrams.tripletId(id0, id1, id2, totalMolecules);
            double betaU = -beta*energy(molecules, r2);
            lastValue[tripletID] = betaU;
            return calcF(betaU);
        }
        double betaUSum = 0;
        // compute sum of triplet energies for all triplets composed of our molecules
        for (int i=0; i<nMolecules-2; i++) {
            int id0 = molecules.getMolecule(i).getIndex();
            for (int j=i+1; j<nMolecules-1; j++) {
                int id1 = molecules.getMolecule(j).getIndex();
                for (int k=j+1; k<nMolecules; k++) {
                    int id2 = molecules.getMolecule(k).getIndex();
                    int tripletID = VirialDiagrams.tripletId(id0, id1, id2, totalMolecules);
                    betaUSum += lastValue[tripletID];
                }
            }
        }
        return calcF(betaUSum);
    }

    /**
     * Returns exp(x)-1 computed directly or using a series expansion
     */
    public double calcF(double x) {
        if (Math.abs(x) < 0.01) {
            return x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        }
        return Math.exp(x)-1;
    }
    
    protected abstract double energy(IMoleculeList molecules, double[] r2);

    public void setBox(IBox box) {
        if (lastValue != null) return;
        IMoleculeList molecules = box.getMoleculeList();
        totalMolecules = molecules.getMoleculeCount();
        if (totalMolecules<3) throw new RuntimeException("need at least 3");
        int nTriplets = VirialDiagrams.tripletId(totalMolecules-3, totalMolecules-2, totalMolecules-1, totalMolecules)+1;
        lastValue = new double[nTriplets];
    }
}
