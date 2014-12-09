/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IMoleculeList;
import etomica.atom.MoleculeArrayList;

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
    protected final MoleculeArrayList molecules;

    public ClusterWheatleyMultibody(int nPoints, MayerFunction f, MayerFunctionNonAdditive fMulti) {
        super(nPoints, f, 1e-12);
        this.fMulti = fMulti;
        moleculeIndices = new int[nPoints];
        r2 = new double[nPoints*(nPoints-1)/2];
        // XXX we don't have a clusterBD for multibody.
        // set it to null so the failure isn't silent
        clusterBD = null;
        molecules = new MoleculeArrayList(nPoints);
    }

    protected void calcFullFQ(BoxCluster box) {
        super.calcFullFQ(box);
        int nf = 1<<n;
        IMoleculeList boxMolecules = box.getMoleculeList();
        // FQ[i] now contains the exp(-bU2) where U2 is the pair-wise energy for set i.
        // we need to go around and add the non-additive energy for each set.
        for (int i=3; i<nf; i++) {
            if (fQ[i] == 0) continue; // pair e-bonds already made this 0
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)) continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            int l = 0;
            molecules.clear();
            for (int a=0; a<n; a++) {
                if ((i & (1<<a)) != 0) {
                    moleculeIndices[l] = a;
                    molecules.add(boxMolecules.getMolecule(a));
                    l++;
                }
            }
            int ll = 0;
            for (int a=0; a<l-1; a++) {
                for (int b=0; b<l; b++) {
                    r2[ll] = box.getCPairSet().getr2(moleculeIndices[a],moleculeIndices[b]);
                    ll++;
                }
            }
            fQ[i] *= fMulti.f(molecules, l, moleculeIndices, r2, beta)+1;
        }
    }
}
