/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.potential.Potential2Soft;
import etomica.space.Space;
import etomica.space.Vector;

public class PotentialMoleculePairPI implements IPotentialMolecular {
    protected final Space space;
    protected Box box;
    protected final Potential2Soft atomPotential;
    protected final double[][] uBeads;
    protected boolean[] pairDone;
    protected long lastPairID, lastBoxIndex;

    public PotentialMoleculePairPI(Space space, Potential2Soft atomPotential, int beadFac) {
        this.space = space;
        this.atomPotential = atomPotential;
        uBeads = new double[beadFac][0];
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0), molecules.get(1));
    }

    public double energy(IMolecule molecule1, IMolecule molecule2) {

        CoordinatePairSet cPairs = ((BoxCluster)box).getCPairSet();
        long thisCPairID = cPairs.getID();
        int beadFac = uBeads.length;
        int mIndex0 = molecule1.getIndex();
        int mIndex1 = molecule2.getIndex();
        int nMolecules = box.getMoleculeList().size();
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
            IAtomList atoms1 = molecule1.getChildList();
            IAtomList atoms2 = molecule2.getChildList();
            for (int i=0; i<atoms1.size(); i++) {
                Vector dr = space.makeVector();
                IAtom atom1 = atoms1.get(i);
                IAtom atom2 = atoms2.get(i);
                dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                uBeads[i%beadFac][thisPairID] += atomPotential.u(dr, atom1, atom2);
            }
            pairDone[thisPairID] = true;
        }
        double uFull = 0;
        for (int i=0; i<uBeads.length; i++) {
            uFull += uBeads[i][thisPairID];
        }
        return uFull;

    }

    @Override
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public void setBox(Box box) {
        this.box = box;
    }

    @Override
    public int nBody() {
        return 2;
    }

    public class PotentialMoleculePISkip implements IPotentialMolecular {

        public PotentialMoleculePISkip(int offset) {
            this.offset = offset;
        }

        public double getRange() {
            return PotentialMoleculePairPI.this.getRange();
        }

        public void setBox(Box box) {
            PotentialMoleculePairPI.this.setBox(box);
        }

        public int nBody() {
            return 2;
        }

        public double energy(IMoleculeList molecules) {
            // poke the potential for the full energy
            // if it has seen these molecules before, it will return without
            // recomputing anything.  if not, it will compute the contributions
            // from all beads
            PotentialMoleculePairPI.this.energy(molecules);
            int mIndex0 = molecules.get(0).getIndex();
            int mIndex1 = molecules.get(1).getIndex();
            int nMolecules = box.getMoleculeList().size();
            int thisPairID = (2*nMolecules-mIndex0-1)*mIndex0/2 + (mIndex1-mIndex0-1);
            // now grab the contribution we care about.
            return PotentialMoleculePairPI.this.uBeads[offset][thisPairID];
        }

        protected final int offset;
    }
}
