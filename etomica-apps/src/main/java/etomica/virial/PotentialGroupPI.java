/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IPotentialMolecular;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialGroup;

/**
 * PotentialGroup that calculates the energy between PI molecule pairs.
 * This class remembers the contribution from each set of beads and can return
 * those energies without additional computation once the full energy has been
 * calculated.
 * 
 * @author Andrew Schultz
 */
public class PotentialGroupPI extends PotentialGroup {

    public PotentialGroupPI(int beadFac) {
        super(2);
        uBeads = new double[beadFac][0];
    }

    public void setBox(Box box) {
        super.setBox(box);
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        int nPairs = nMolecules*(nMolecules-1)/2;
        if (pairDone != null && pairDone.length == nPairs) return;
        pairDone = new boolean[nPairs];
        for (int i=0; i<uBeads.length; i++) {
            uBeads[i] = new double[nPairs];
        }
    }
    
    public double energy(IMoleculeList molecules) {
        CoordinatePairSet cPairs = ((BoxCluster)box).getCPairSet();
        long thisCPairID = cPairs.getID();
        int beadFac = uBeads.length;
        int mIndex0 = molecules.getMolecule(0).getIndex();
        int mIndex1 = molecules.getMolecule(1).getIndex();
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        int thisPairID = (2*nMolecules-mIndex0-1)*mIndex0/2 + (mIndex1-mIndex0-1);
        if (thisCPairID != lastPairID || box.getIndex() != lastBoxIndex) {
            for (int i=0; i<pairDone.length; i++) {
                pairDone[i] = false;
            }
        }
        if (thisCPairID != lastPairID || box.getIndex() != lastBoxIndex || !pairDone[thisPairID]) {
            // we haven't seen this pair of molecules before
            // recalculate everything, remembering the contributions for 0,2,4,6...
            // and 1,3,5,7... separately.
            lastPairID = thisCPairID;
            lastBoxIndex = box.getIndex();
            for (int j=0; j<beadFac; j++) {
                uBeads[j][thisPairID] = 0;
            }
            for (PotentialLinker link=first; link!= null; link=link.next) { 
                if(!link.enabled) continue;
                link.iterator.setBasis(molecules);
                link.iterator.reset();
                int i=0;
                for (IAtomList atoms = link.iterator.next(); atoms != null; atoms = link.iterator.next()) {
                    uBeads[i%beadFac][thisPairID] += link.potential.energy(atoms);
                    i++;
                }
            }
            pairDone[thisPairID] = true;
        }
        double uFull = 0;
        for (int i=0; i<uBeads.length; i++) {
            uFull += uBeads[i][thisPairID];
        }
        return uFull;
    }

    public class PotentialGroupPISkip implements IPotentialMolecular {

        public PotentialGroupPISkip(int offset) {
            this.offset = offset;
        }
        
        public double getRange() {
            return PotentialGroupPI.this.getRange();
        }

        public void setBox(Box box) {
            PotentialGroupPI.this.setBox(box);
        }

        public int nBody() {
            return 2;
        }

        public double energy(IMoleculeList molecules) {
            // poke the potential for the full energy
            // if it has seen these molecules before, it will return without
            // recomputing anything.  if not, it will compute the contributions
            // from all beads
            PotentialGroupPI.this.energy(molecules);
            int mIndex0 = molecules.getMolecule(0).getIndex();
            int mIndex1 = molecules.getMolecule(1).getIndex();
            int nMolecules = box.getMoleculeList().getMoleculeCount();
            int thisPairID = (2*nMolecules-mIndex0-1)*mIndex0/2 + (mIndex1-mIndex0-1);
            // now grab the contribution we care about.
            return PotentialGroupPI.this.uBeads[offset][thisPairID];
        }
        
        protected final int offset;
    }
    
    protected final double[][] uBeads;
    protected boolean[] pairDone;
    protected long lastPairID, lastBoxIndex;
}
