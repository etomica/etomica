/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialHard;
import etomica.potential.Potential2Soft;
import etomica.space.Vector;
import etomica.util.collections.Int2IntHash;

/**
 * Simple (all atoms are neighbors of all atoms) neighbor manager that keeps
 * track of state for pairs of atoms.  We use an int-int hashmap, which isn't
 * super-fast, but should be fine and faster than the alternatives.
 */
public class NeighborManagerSimpleHard extends NeighborManagerSimple implements NeighborManagerHard {

    public Int2IntHash[] stateHash;
    protected IPotentialAtomic[][] pairPotentials;

    public NeighborManagerSimpleHard(Box box) {
        super(box);
    }

    @Override
    public NeighborIteratorHard makeNeighborIterator() {
        return new NeighborIteratorSimpleHard(this);
    }

    @Override
    public void setPairPotentials(Potential2Soft[][] potentials) {
        pairPotentials = potentials;
    }

    @Override
    public void init() {
        super.init();
        reset();
    }

    protected void reset() {
        IAtomList atoms = box.getLeafList();
        stateHash = new Int2IntHash[atoms.size()];
        for (int i = 0; i < stateHash.length; i++) {
            stateHash[i] = new Int2IntHash(10);
            IAtomKinetic iAtom = (IAtomKinetic) atoms.get(i);
            Vector ri = iAtom.getPosition();
            IPotentialAtomic[] iPotentials = pairPotentials[iAtom.getType().getIndex()];
            for (int j = i + 1; j < stateHash.length; j++) {
                IAtomKinetic jAtom = (IAtomKinetic) atoms.get(j);
                IPotentialHard p2 = (IPotentialHard) iPotentials[jAtom.getType().getIndex()];
                if (p2 == null) continue;
                Vector rj = jAtom.getPosition();
                Vector dr = Vector.d(ri.getD());
                dr.Ev1Mv2(rj, ri);
                box.getBoundary().nearestImage(dr);
                int state = p2.getState(iAtom, jAtom, dr);
                if (state >= 0) stateHash[i].put(j, state);
            }
        }
    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return new IntegratorListener() {
            @Override
            public void integratorInitialized(IntegratorEvent e) {
                init();
            }

            @Override
            public void integratorStepStarted(IntegratorEvent e) {
            }

            @Override
            public void integratorStepFinished(IntegratorEvent e) {
            }
        };
    }

    @Override
    public void setPairState(int i, int j, int state) {
        Int2IntHash iHash = stateHash[i];
        int oldState = iHash.get(j);
        if (state == oldState || (state < 0 && oldState < 0)) return;
        if (state < 0) {
            iHash.remove(j);
        } else {
            iHash.put(j, state);
        }
    }

}
