/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Vector;

class NeighborIteratorSimpleHard implements NeighborIteratorHard {

    private final NeighborManagerSimpleHard neighborManager;

    public NeighborIteratorSimpleHard(NeighborManagerSimpleHard neighborManager) {
        this.neighborManager = neighborManager;
    }

    @Override
    public void iterUpNeighbors(int i, NeighborConsumerHard consumer, double falseTime) {
        IAtomList atoms = neighborManager.box.getLeafList();
        Vector rij = neighborManager.box.getSpace().makeVector();
        Vector ri = atoms.get(i).getPosition();
        Vector vi = ((IAtomKinetic) atoms.get(i)).getVelocity();
        Vector vij = neighborManager.box.getSpace().makeVector();
        for (int j = i + 1; j < atoms.size(); j++) {
            if (neighborManager.pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                continue;
            rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
            vij.Ev1Mv2(((IAtomKinetic) atoms.get(j)).getVelocity(), vi);
            rij.PEa1Tv1(falseTime, vij);
            neighborManager.box.getBoundary().nearestImage(rij);
            consumer.acceptHard(atoms.get(j), rij, neighborManager.stateHash[i].get(j));
        }
    }

    @Override
    public void iterUpNeighbors(int i, NeighborConsumer consumer) {
        IAtomList atoms = neighborManager.box.getLeafList();
        Vector rij = neighborManager.box.getSpace().makeVector();
        Vector ri = atoms.get(i).getPosition();
        for (int j = i + 1; j < atoms.size(); j++) {
            if (neighborManager.pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                continue;
            rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
            neighborManager.box.getBoundary().nearestImage(rij);
            consumer.accept(atoms.get(j), rij, 0);
        }
    }

    @Override
    public void iterDownNeighbors(int i, NeighborConsumerHard consumer, double falseTime) {
        IAtomList atoms = neighborManager.box.getLeafList();
        Vector rij = neighborManager.box.getSpace().makeVector();
        Vector ri = atoms.get(i).getPosition();
        Vector vi = ((IAtomKinetic) atoms.get(i)).getVelocity();
        Vector vij = neighborManager.box.getSpace().makeVector();
        for (int j = 0; j < i; j++) {
            if (neighborManager.pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                continue;
            rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
            vij.Ev1Mv2(((IAtomKinetic) atoms.get(j)).getVelocity(), vi);
            rij.PEa1Tv1(falseTime, vij);
            neighborManager.box.getBoundary().nearestImage(rij);
            consumer.acceptHard(atoms.get(j), rij, neighborManager.stateHash[j].get(i));
        }
    }

    @Override
    public void iterDownNeighbors(int i, NeighborConsumer consumer) {
        IAtomList atoms = neighborManager.box.getLeafList();
        Vector rij = neighborManager.box.getSpace().makeVector();
        Vector ri = atoms.get(i).getPosition();
        for (int j = 0; j < i; j++) {
            if (neighborManager.pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                continue;
            rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
            neighborManager.box.getBoundary().nearestImage(rij);
            consumer.accept(atoms.get(j), rij, 0);
        }
    }

    @Override
    public void iterAllNeighbors(int i, NeighborConsumerHard consumer, double falseTime) {
        IAtomList atoms = neighborManager.box.getLeafList();
        Vector rij = neighborManager.box.getSpace().makeVector();
        Vector ri = atoms.get(i).getPosition();
        Vector vi = ((IAtomKinetic) atoms.get(i)).getVelocity();
        Vector vij = neighborManager.box.getSpace().makeVector();
        for (int j = 0; j < atoms.size(); j++) {
            if (j == i) continue;
            if (neighborManager.pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                continue;
            rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
            vij.Ev1Mv2(((IAtomKinetic) atoms.get(j)).getVelocity(), vi);
            rij.PEa1Tv1(falseTime, vij);
            neighborManager.box.getBoundary().nearestImage(rij);
            consumer.acceptHard(atoms.get(j), rij, neighborManager.stateHash[i].get(j));
        }
    }

    @Override
    public void iterAllNeighbors(int i, NeighborConsumer consumer) {
        IAtomList atoms = neighborManager.box.getLeafList();
        Vector rij = neighborManager.box.getSpace().makeVector();
        Vector ri = atoms.get(i).getPosition();
        for (int j = 0; j < atoms.size(); j++) {
            if (j == i) continue;
            if (neighborManager.pairPotentials[atoms.get(i).getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                continue;
            rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
            neighborManager.box.getBoundary().nearestImage(rij);
            consumer.accept(atoms.get(j), rij, 0);
        }
    }

    @Override
    public double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {
        IAtomList atoms = neighborManager.box.getLeafList();
        Vector rij = neighborManager.box.getSpace().makeVector();
        Vector ri = atom1.getPosition();
        double sum = 0;
        for (int j = 0; j < atoms.size(); j++) {
            if (j == atom1.getLeafIndex()) continue;
            if (neighborManager.pairPotentials[atom1.getType().getIndex()][atoms.get(j).getType().getIndex()] == null)
                continue;
            rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
            neighborManager.box.getBoundary().nearestImage(rij);
            sum += consumer.accept(atom1, atoms.get(j), rij, 0);
        }
        return sum;
    }
}
