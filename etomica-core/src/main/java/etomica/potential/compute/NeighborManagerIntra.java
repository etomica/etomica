/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.potential.BondingInfo;
import etomica.space.Vector;

/**
 * NeighborManager that brute-force loops over all appropriate pairs.
 */
public class NeighborManagerIntra implements NeighborManager {

    protected final Box box;
    protected final BondingInfo bondingInfo;

    public NeighborManagerIntra(Box box, BondingInfo bondingInfo) {
        this.box = box;
        this.bondingInfo = bondingInfo;
    }

    public Box getBox() {
        return box;
    }

    @Override
    public BondingInfo getBondingInfo() {
        return bondingInfo;
    }

    @Override
    public NeighborIterator makeNeighborIterator() {
        return new NeighborIterator() {
            @Override
            public void iterUpNeighbors(int i, NeighborConsumer consumer) {
                IAtomList atoms = box.getLeafList();
                IAtom iAtom = atoms.get(i);
                IMolecule iMolecule = iAtom.getParentGroup();
                Vector rij = box.getSpace().makeVector();
                Vector ri = iAtom.getPosition();
                IAtomList childAtoms = iMolecule.getChildList();
                for (int j = iAtom.getIndex()+1; j<childAtoms.size(); j++) {
                    IAtom jAtom = childAtoms.get(j);
                    if (bondingInfo.skipBondedPair(false, iAtom, jAtom)) continue;
                    rij.Ev1Mv2(jAtom.getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(jAtom, rij);
                }
            }

            @Override
            public void iterDownNeighbors(int i, NeighborConsumer consumer) {
                IAtomList atoms = box.getLeafList();
                IAtom iAtom = atoms.get(i);
                IMolecule iMolecule = iAtom.getParentGroup();
                Vector rij = box.getSpace().makeVector();
                Vector ri = iAtom.getPosition();
                IAtomList childAtoms = iMolecule.getChildList();
                for (int j = 0; j<iAtom.getIndex(); j++) {
                    IAtom jAtom = childAtoms.get(j);
                    if (jAtom.getParentGroup() != iMolecule) continue;
                    if (bondingInfo.skipBondedPair(false, iAtom, jAtom)) continue;
                    rij.Ev1Mv2(jAtom.getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(jAtom, rij);
                }
            }

            @Override
            public void iterAllNeighbors(int i, NeighborConsumer consumer) {
                IAtomList atoms = box.getLeafList();
                IAtom iAtom = atoms.get(i);
                IMolecule iMolecule = iAtom.getParentGroup();
                Vector rij = box.getSpace().makeVector();
                Vector ri = iAtom.getPosition();
                IAtomList childAtoms = iMolecule.getChildList();
                for (int j = 0; j<childAtoms.size(); j++) {
                    if (i==j) continue;
                    IAtom jAtom = childAtoms.get(j);
                    if (jAtom.getParentGroup() != iMolecule) continue;
                    if (bondingInfo.skipBondedPair(false, iAtom, jAtom)) continue;
                    rij.Ev1Mv2(jAtom.getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(jAtom, rij);
                }
            }

            @Override
            public double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {
                IAtomList atoms = box.getLeafList();
                IMolecule iMolecule = atom1.getParentGroup();
                Vector rij = box.getSpace().makeVector();
                Vector ri = atom1.getPosition();
                double sum = 0;
                IAtomList childAtoms = iMolecule.getChildList();
                for (int j = 0; j<childAtoms.size(); j++) {
                    if (atom1.getIndex()==j) continue;
                    IAtom jAtom = childAtoms.get(j);
                    if (bondingInfo.skipBondedPair(false, atom1, jAtom)) continue;
                    rij.Ev1Mv2(jAtom.getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    sum += consumer.accept(atom1, jAtom, rij);
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
        return new IntegratorListener() {};
    }

    @Override
    public void updateAtom(IAtom atom) {

    }
}
