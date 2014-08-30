package etomica.virial;

/**
 * This class uses Wheatley's recursion approach to calculating all biconnected
 * diagrams, but adds in non-additive contributions.
 * 
 * @author Andrew Schultz
 */
public class ClusterWheatleyMultibody extends ClusterWheatleySoft {

    protected final MayerFunctionNonAdditive fMulti;
    protected final int[] moleculeIndices;
    protected final double[] r2;

    public ClusterWheatleyMultibody(int nPoints, MayerFunction f, MayerFunctionNonAdditive fMulti) {
        super(nPoints, f, 1e-12);
        this.fMulti = fMulti;
        moleculeIndices = new int[nPoints];
        r2 = new double[nPoints*(nPoints-1)/2];
        // clusterBD is going to fail completely here... set its temperature to
        // 0 so the failure isn't silent 
        clusterBD.setTemperature(0);
    }

    protected void calcFullFQ(BoxCluster box) {
        super.calcFullFQ(box);
        int nf = 1<<n;
        // FQ[i] now contains the exp(-bU2) where U2 is the pair-wise energy for set i.
        // we need to go around and add the non-additive energy for each set.
        for (int i=3; i<nf; i++) {
            if (fQ[i] == 0) continue; // pair e-bonds already made this 0
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)) continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            int l = 0;
            for (int a=0; a<n; a++) {
                if ((i & (1<<a)) != 0) {
                    moleculeIndices[l] = a;
                    l++;
                }
            }
            int ll = 0;
            for (int a=0; a<l-1; a++) {
                for (int b=0; b<l; b++) {
                    r2[ll] = box.getCPairSet().getr2(a,b);
                }
            }
            fQ[i] *= fMulti.f(null, l, moleculeIndices, r2, beta);
        }
    }
}
