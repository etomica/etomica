/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.virial.cluster.VirialDiagrams;


public class PotentialCommonAtomic implements IPotentialAtomic {
    protected IPotentialAtomic pot;
    public PotentialCommonAtomic(IPotentialAtomic pot1) {
        this.pot = pot1;
        nBody = pot.nBody();
    }

    public double getRange() {
        return pot.getRange();
    }

    public void setBox(Box box) {
        this.box1 = box;
        pot.setBox(box);
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        int nPairs = nMolecules*(nMolecules-1)/2;
        if (pairDone != null && pairDone.length == nPairs && tripletDone != null) return;
        pairDone = new boolean[nPairs];
        uPair = new double[nPairs];
        if (nBody == 3) {
            int nTriplets = VirialDiagrams.tripletId(nMolecules-3, nMolecules-2, nMolecules-1, nMolecules)+1;
            tripletDone = new boolean[nTriplets];
            uTriplet = new double[nTriplets];
        }
    }

    public int nBody() {
        return nBody;
    }

    public double energy(IAtomList atoms) {
        CoordinatePairSet cPairs = ((BoxCluster)box1).getCPairSet();
        long thisCPairID = cPairs.getID();        
        int aIndex0 = atoms.get(0).getLeafIndex();
        int aIndex1 = atoms.get(1).getLeafIndex();
        int nMolecules = box1.getLeafList().size();
        int thisPairID = (2*nMolecules-aIndex0-1)*aIndex0/2 + (aIndex1-aIndex0-1);
        if (thisCPairID != lastPairID || box1.getIndex() != lastBoxIndex) {
            if (nBody == 2) {
                for (int i=0; i<pairDone.length; i++) {
                    pairDone[i] = false;
                }
            }
            if (nBody == 3) {
                for (int i=0; i<tripletDone.length; i++) {
                    tripletDone[i] = false;
                }
            }
        }
        if (nBody == 2) {            
            if (!pairDone[thisPairID]) {
                // we haven't seen this pair of molecules before
                lastPairID = thisCPairID;
                lastBoxIndex = box1.getIndex();
                uPair[thisPairID] = pot.energy(atoms);
                pairDone[thisPairID] = true;                
            }
            return uPair[thisPairID];
        }
        else if (nBody == 3) {
            int aIndex2 = atoms.get(2).getLeafIndex();
            int tripletID = VirialDiagrams.tripletId(aIndex0, aIndex1, aIndex2, nMolecules);            
            if (!tripletDone[tripletID]) {
                // we haven't seen this pair of molecules before
                lastPairID = thisCPairID;
                lastBoxIndex = box1.getIndex();
                uTriplet[tripletID] = pot.energy(atoms);
                tripletDone[tripletID] = true;
            }
            return uTriplet[tripletID];
        }
        return 0;
    }
    protected boolean[] pairDone, tripletDone;
    protected long lastPairID, lastBoxIndex;
    protected Box box1;
    protected double[] uPair, uTriplet;
    protected int nBody;
}
