/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial;

import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.potential.Potential3Soft;
import etomica.space.Space;
import etomica.virial.cluster.VirialDiagrams;

public class PotentialMolecule3PI implements IPotentialMolecular {
    protected final Space space;
    protected final int nMolecules;
    protected final Potential3Soft atomPotential;
    protected final double[][] uBeads;

    public PotentialMolecule3PI(Space space, Potential3Soft atomPotential, int beadFac, int nMolecules) {
        this.space = space;
        this.nMolecules = nMolecules;
        this.atomPotential = atomPotential;
        uBeads = new double[beadFac][0];
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0), molecules.get(1), molecules.get(2));
    }

    public double energy(IMolecule molecule1, IMolecule molecule2, IMolecule molecule3) {

        int beadFac = uBeads.length;

        int id1 = molecule1.getIndex();
        int id2 = molecule2.getIndex();
        int id3 = molecule3.getIndex();
        int tripletID = VirialDiagrams.tripletId(id1, id2, id3, nMolecules);

        // TODO do we get invoked for the same pair without anything moving?
        // TODO can we avoid recomputing everything every time?
        // recalculate everything, remembering the contributions for 0,2,4,6...
        // and 1,3,5,7... separately.
        for (int j=0; j<beadFac; j++) {
            uBeads[j][tripletID] = 0;
        }
        IAtomList atoms1 = molecule1.getChildList();
        IAtomList atoms2 = molecule2.getChildList();
        IAtomList atoms3 = molecule3.getChildList();
        for (int i=0; i<atoms1.size(); i++) {
            double r212 = atoms1.get(i).getPosition().Mv1Squared(atoms2.get(i).getPosition());
            double r213 = atoms1.get(i).getPosition().Mv1Squared(atoms3.get(i).getPosition());
            double r223 = atoms2.get(i).getPosition().Mv1Squared(atoms3.get(i).getPosition());
            uBeads[i%beadFac][tripletID] += atomPotential.u(r212, r213, r223);
        }
        double uFull = 0;
        for (int i=0; i<uBeads.length; i++) {
            uFull += uBeads[i][tripletID];
        }
        return uFull;

    }

    @Override
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public class PotentialMolecule3PISkip implements IPotentialMolecular {

        public PotentialMolecule3PISkip(int offset) {
            this.offset = offset;
        }

        public double getRange() {
            return PotentialMolecule3PI.this.getRange();
        }

        public double energy(IMoleculeList molecules) {
            // poke the potential for the full energy
            // if it has seen these molecules before, it will return without
            // recomputing anything.  if not, it will compute the contributions
            // from all beads
            PotentialMolecule3PI.this.energy(molecules);
            int id0 = molecules.get(0).getIndex();
            int id1 = molecules.get(1).getIndex();
            int id2 = molecules.get(2).getIndex();
            int tripletID = VirialDiagrams.tripletId(id0, id1, id2, nMolecules);

            return PotentialMolecule3PI.this.uBeads[offset][tripletID];
        }

        protected final int offset;
    }
}
