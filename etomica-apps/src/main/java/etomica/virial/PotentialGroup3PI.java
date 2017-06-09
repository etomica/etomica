/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.api.IPotentialMolecular;
import etomica.potential.PotentialGroup;
import etomica.virial.cluster.VirialDiagrams;

/**
 * PotentialGroup that calculates the energy between PI molecules triplets.
 * This class remembers the contribution from each set of beads and can return
 * those energies without additional computation once the full enegy has been
 * calculated.
 * 
 * @see PotentialGroupPI
 * @author Andrew Schultz
 */
public class PotentialGroup3PI extends PotentialGroup {

    public PotentialGroup3PI(int beadFac) {
        super(3);
        uBeads = new double[beadFac][0];
    }

    public void setBox(Box box) {
        super.setBox(box);
        if (tripletDone != null) return;
        int nMolecules = box.getMoleculeList().getMoleculeCount();

        int nTriplets = VirialDiagrams.tripletId(nMolecules-3, nMolecules-2, nMolecules-1, nMolecules)+1;

        tripletDone = new boolean[nTriplets];
        for (int i=0; i<uBeads.length; i++) {
            uBeads[i] = new double[nTriplets];
        }
    }
    
    public double energy(IMoleculeList molecules) {
        CoordinatePairSet cPairs = ((BoxCluster)box).getCPairSet();
        long thisCPairID = cPairs.getID();
        int beadFac = uBeads.length;
        int nMolecules = box.getMoleculeList().getMoleculeCount();

        if (thisCPairID != lastPairID || box.getIndex() != lastBoxIndex) {
            for (int i=0; i<tripletDone.length; i++) {
                tripletDone[i] = false;
            }
        }

        int id0 = molecules.getMolecule(0).getIndex();
        int id1 = molecules.getMolecule(1).getIndex();
        int id2 = molecules.getMolecule(2).getIndex();
        int tripletID = VirialDiagrams.tripletId(id0, id1, id2, nMolecules);

        if (thisCPairID != lastPairID || box.getIndex() != lastBoxIndex || !tripletDone[tripletID]) {
            lastPairID = thisCPairID;
            lastBoxIndex = box.getIndex();
            for (int j=0; j<beadFac; j++) {
                uBeads[j][tripletID] = 0;
            }
            for (PotentialLinker link=first; link!= null; link=link.next) { 
                if(!link.enabled) continue;
                link.iterator.setBasis(molecules);
                link.iterator.reset();
                int i=0;
                for (IAtomList atoms = link.iterator.next(); atoms != null; atoms = link.iterator.next()) {
                    uBeads[i%beadFac][tripletID] += link.potential.energy(atoms);
                    i++;
                }
            }
            tripletDone[tripletID] = true;
        }
        double uFull = 0;
        for (int i=0; i<uBeads.length; i++) {
            uFull += uBeads[i][tripletID];
        }
        return uFull;
    }

    public class PotentialGroup3PISkip implements IPotentialMolecular {

        public PotentialGroup3PISkip(int offset) {
            this.offset = offset;
        }
        
        public double getRange() {
            return PotentialGroup3PI.this.getRange();
        }

        public void setBox(Box box) {
            PotentialGroup3PI.this.setBox(box);
        }

        public int nBody() {
            return 2;
        }

        public double energy(IMoleculeList molecules) {
            PotentialGroup3PI.this.energy(molecules);

            int nMolecules = box.getMoleculeList().getMoleculeCount();
            int id0 = molecules.getMolecule(0).getIndex();
            int id1 = molecules.getMolecule(1).getIndex();
            int id2 = molecules.getMolecule(2).getIndex();
            int tripletID = VirialDiagrams.tripletId(id0, id1, id2, nMolecules);

            return PotentialGroup3PI.this.uBeads[offset][tripletID];
        }
        
        protected final int offset;
    }
    
    protected final double[][] uBeads;
    protected boolean[] tripletDone;
    protected long lastPairID, lastBoxIndex;
}
