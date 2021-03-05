/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialHard;
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
    public NeighborIterator makeNeighborIterator() {
        return new NeighborIterator() {
            @Override
            public void iterUpNeighbors(int i, NeighborConsumer consumer) {
                final NeighborConsumerHard consumerHard = consumer instanceof NeighborConsumerHard ? (NeighborConsumerHard) consumer : null;
                IAtomList atoms = box.getLeafList();
                Vector rij = box.getSpace().makeVector();
                Vector ri = atoms.get(i).getPosition();
                for (int j = i + 1; j < atoms.size(); j++) {
                    if (pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                        continue;
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    if (consumerHard == null) {
                        // we're like just computing energy
                        consumer.accept(atoms.get(j), rij);
                    } else {
                        consumerHard.accept(atoms.get(j), rij, stateHash[i].get(j));
                    }
                }
            }

            @Override
            public void iterDownNeighbors(int i, NeighborConsumer consumer) {
                final NeighborConsumerHard consumerHard = consumer instanceof NeighborConsumerHard ? (NeighborConsumerHard) consumer : null;
                IAtomList atoms = box.getLeafList();
                Vector rij = box.getSpace().makeVector();
                Vector ri = atoms.get(i).getPosition();
                for (int j = 0; j < i; j++) {
                    if (pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                        continue;
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    if (consumerHard == null) {
                        // we're like just computing energy
                        consumer.accept(atoms.get(j), rij);
                    } else {
                        consumerHard.accept(atoms.get(j), rij, stateHash[j].get(i));
                    }
                }
            }

            @Override
            public void iterAllNeighbors(int i, NeighborConsumer consumer) {
                final NeighborConsumerHard consumerHard = consumer instanceof NeighborConsumerHard ? (NeighborConsumerHard) consumer : null;
                IAtomList atoms = box.getLeafList();
                Vector rij = box.getSpace().makeVector();
                Vector ri = atoms.get(i).getPosition();
                for (int j = 0; j < atoms.size(); j++) {
                    if (j == i) continue;
                    if (pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                        continue;
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(atoms.get(j), rij);
                    if (consumerHard == null) {
                        // we're like just computing energy
                        consumer.accept(atoms.get(j), rij);
                    } else {
                        consumerHard.accept(atoms.get(j), rij, stateHash[i].get(j));
                    }
                }
            }

            @Override
            public double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {
                IAtomList atoms = box.getLeafList();
                Vector rij = box.getSpace().makeVector();
                Vector ri = atom1.getPosition();
                double sum = 0;
                for (int j = 0; j < atoms.size(); j++) {
                    if (j == atom1.getLeafIndex()) continue;
                    if (pairPotentials[atom1.getType().getIndex()][atoms.get(j).getType().getIndex()] == null) continue;
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    sum += consumer.accept(atom1, atoms.get(j), rij);
                }
                return sum;
            }
        };
    }

    @Override
    public void setPairPotentials(IPotentialAtomic[][] potentials) {
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
