/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.osmoticvirial;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.nbr.cell.NeighborCellManagerFasterer;
import etomica.potential.BondingInfo;
import etomica.potential.IPotentialAtomic;
import etomica.potential.compute.NeighborIterator;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

public class NeighborManagerCellMixed extends NeighborCellManagerFasterer {

    protected final AtomType cellType;
    protected IPotentialAtomic[][] pairPotentials;

    public NeighborManagerCellMixed(SpeciesManager sm, Box box, int cellRange, BondingInfo bondingInfo, AtomType cellType) {
        super(sm, box, cellRange, bondingInfo);
        this.cellType = cellType;
    }

    @Override
    public void setPairPotentials(IPotentialAtomic[][] potentials) {
        pairPotentials = potentials;
        setPotentialRange(0);
    }

    @Override
    public void setPotentialRange(double newRange) {
        if (pairPotentials == null) return;
        IPotentialAtomic p = pairPotentials[cellType.getIndex()][cellType.getIndex()];
        if (p != null) {
            super.setPotentialRange(p.getRange());
        }
    }

    public NeighborIterator makeNeighborIterator() {
        NeighborIterator nic = super.makeNeighborIterator();
        IAtomList atoms = box.getLeafList();
        return new NeighborIterator() {
            @Override
            public void iterUpNeighbors(int iAtom, NeighborConsumer consumer) {
                boolean isCellType = box.getLeafList().get(iAtom).getType() == cellType;
                if (isCellType) {
                    nic.iterUpNeighbors(iAtom, new NeighborConsumer() {
                        @Override
                        public void accept(IAtom jAtom, Vector rij) {
                            if (jAtom.getType() == cellType) consumer.accept(jAtom, rij);
                        }
                    });
                }
                Vector ri = atoms.get(iAtom).getPosition();
                for (int j = iAtom+1; j<atoms.size(); j++) {
                    if (isCellType && atoms.get(j).getType() == cellType) continue;
                    Vector rij = box.getSpace().makeVector();
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(atoms.get(j), rij);
                }
            }

            @Override
            public void iterDownNeighbors(int iAtom, NeighborConsumer consumer) {
                boolean isCellType = box.getLeafList().get(iAtom).getType() == cellType;
                if (isCellType) {
                    nic.iterDownNeighbors(iAtom, new NeighborConsumer() {
                        @Override
                        public void accept(IAtom jAtom, Vector rij) {
                            if (jAtom.getType() == cellType) consumer.accept(jAtom, rij);
                        }
                    });
                }
                Vector ri = atoms.get(iAtom).getPosition();
                for (int j = 0; j<iAtom; j++) {
                    if (isCellType && atoms.get(j).getType() == cellType) continue;
                    Vector rij = box.getSpace().makeVector();
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(atoms.get(j), rij);
                }
            }

            @Override
            public void iterAllNeighbors(int iAtom, NeighborConsumer consumer) {
                boolean isCellType = box.getLeafList().get(iAtom).getType() == cellType;
                if (isCellType) {
                    nic.iterAllNeighbors(iAtom, new NeighborConsumer() {
                        @Override
                        public void accept(IAtom jAtom, Vector rij) {
                            if (jAtom.getType() == cellType) consumer.accept(jAtom, rij);
                        }
                    });
                }

                Vector ri = atoms.get(iAtom).getPosition();
                for (int j = 0; j<atoms.size(); j++) {
                    if (j == iAtom) continue;
                    if (isCellType && atoms.get(j).getType() == cellType) continue;
                    Vector rij = box.getSpace().makeVector();
                    rij.Ev1Mv2(atoms.get(j).getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    consumer.accept(atoms.get(j), rij);
                }
            }

            @Override
            public double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {
                double sum = 0;
                boolean isCellType = atom1.getType() == cellType;
                if (isCellType) {
                    sum = nic.iterAndSumAllNeighbors(atom1, new SuperNbrConsumer() {
                        @Override
                        public double accept(IAtom atom1, IAtom jAtom, Vector rij) {
                            if (jAtom.getType() != cellType) return 0;
                            return consumer.accept(atom1, jAtom, rij);
                        }
                    });
                }

                Vector ri = atom1.getPosition();
                for (int j = 0; j<atoms.size(); j++) {
                    if (j == atom1.getLeafIndex()) continue;
                    IAtom jAtom = atoms.get(j);
                    if (isCellType && jAtom.getType() == cellType) continue;
                    Vector rij = box.getSpace().makeVector();
                    rij.Ev1Mv2(jAtom.getPosition(), ri);
                    box.getBoundary().nearestImage(rij);
                    sum += consumer.accept(atom1, jAtom, rij);
                }
                return sum;
            }
        };
    }
}
