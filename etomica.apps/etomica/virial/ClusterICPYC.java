package etomica.virial;


/**
 * Cluster class that computes the correction to the ICPY (incrementally
 * corrected Percus-Yevick) formulation.
 *
 * @author Andrew Schultz
 */
public class ClusterICPYC extends ClusterWheatley {

    protected final boolean allPermutations = true;
    protected int num, tot;
    
    public ClusterICPYC(int nPoints, MayerFunction f) {
        super(nPoints, f);
    }

    protected void calcValue(BoxCluster box) {
        /*
         * The recipe for the incrementally-corrected PY correction is
         * 1. take all biconnected f-diagrams that do not have a bond between 0-1
         * 2. add an e-bond between 0-1
         * 
         * We can get the value for #1 by setting f12=0 and computing all biconnected
         * diagrams.  Then, we multiply by e12.
         */
        double icpycValue = 0;
        tot++;
        
        if (allPermutations && doBiconCheck) {
            // if the full graph is not biconnected, we won't be able to remove a bond and make it biconnected!
            int bit=0;
            for (int i=0; i<n; i++) {
                outDegree[i] = 0;
            }

            for (int i=0; i<n-1; i++) {
                for (int j=i+1; j<n; j++) {
                    boolean fBond = (fQ[(1<<i)|(1<<j)] == 0); 
                    if (fBond) {
                        bondMap.setBit(bit);
                        outDegree[i]++;
                        outDegree[j]++;
                    }
                    else {
                        bondMap.clearBit(bit);
                    }
                    bit++;
                }
            }
            for (int i=0; i<n; i++) {
                if (outDegree[i] < 2) {
                    value = 0;
                    return;
                }
            }
            if (!isBi.check(myGraph)) {
                value = 0;
                return;
            }
        }
        
        for (int i=0; i<n-1; i++) {
            for (int j=i+1; j<n; j++) {
                int k = (1<<j)|(1<<i);
                double ek = fQ[k];
                if (ek == 0) {
                    if (!allPermutations) {
                        value = 0;
                        return;
                    }
                    continue;
                }
                fQ[k] = 1;
                super.calcValue(box);
                fQ[k] = ek;
                icpycValue += value*ek;
                if (!allPermutations) {
                    value = icpycValue;
                    return;
                }
            }
        }
        // we've computed the ICPYC value for each pair of root points.
        // now divide by the number of pairs.
        // this only helps precision a bit
        value = icpycValue / (n*(n-1)/2);
    }
}
