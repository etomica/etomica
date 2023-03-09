/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.BondingInfo;
import etomica.space.Vector;

/**
 * NeighborManager that brute-force loops over all appropriate pairs.
 */
public class NeighborManagerSimple implements NeighborManager {

    protected final Box box;

    public NeighborManagerSimple(Box box) {
        this.box = box;
    }

    public Box getBox() {
        return box;
    }

    @Override
    // this whole class assumes no bonding
    public BondingInfo getBondingInfo() {
        return BondingInfo.noBonding();
    }

    @Override
    public NeighborIterator makeNeighborIterator() {
        return new NeighborIterator() {
            @Override
            public void iterUpNeighbors(int i, NeighborConsumer consumer) {
                IAtomList atoms = box.getLeafList();
                Vector rij = box.getSpace().makeVector();
                Vector ri = atoms.get(i).getPosition();
                for (int j = i + 1; j < atoms.size(); j++) {
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(atoms.get(j), rij, 0);
                }
            }

            @Override
            public void iterDownNeighbors(int i, NeighborConsumer consumer) {
                IAtomList atoms = box.getLeafList();
                Vector rij = box.getSpace().makeVector();
                Vector ri = atoms.get(i).getPosition();
                for (int j = 0; j < i; j++) {
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(atoms.get(j), rij, 0);
                }
            }

            @Override
            public void iterAllNeighbors(int i, NeighborConsumer consumer) {
                IAtomList atoms = box.getLeafList();
                Vector rij = box.getSpace().makeVector();
                Vector ri = atoms.get(i).getPosition();
                for (int j = 0; j < atoms.size(); j++) {
                    if (j == i) continue;
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(atoms.get(j), rij, 0);
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
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    sum += consumer.accept(atom1, atoms.get(j), rij, 0);
                }
                return sum;
            }
        };
    }

    @Override
    public void init() {

    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return new IntegratorListener() {
            @Override
            public void integratorInitialized(IntegratorEvent e) {

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
    public void updateAtom(IAtom atom) {

    }
}
