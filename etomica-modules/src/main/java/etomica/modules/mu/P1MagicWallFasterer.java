/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.potential.IPotentialPair;
import etomica.potential.P1HardFieldGeneric;
import etomica.potential.P2HardGeneric;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManagerHard;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.util.collections.Int2IntHash;

/**
 *
 */

public class P1MagicWallFasterer extends P1HardFieldGeneric {

    protected final NeighborManagerHard neighborManager;
    protected final IPotentialPair[][] potentials;
    protected final Int2IntHash stateHash;

    public P1MagicWallFasterer(double L, PotentialComputePairGeneral potentialMaster, NeighborManagerHard neighborManager) {
        super(0, new double[]{-L / 2, 0, +L / 2}, new double[]{0, 0});
        this.neighborManager = neighborManager;
        this.potentials = potentialMaster.getPairPotentials();
        stateHash = new Int2IntHash(10);
    }

    public double energy(IAtomList a) {
        double e = 0.0;
        // this is probably wrong, unfortunately
        return e;
    }

    @Override
    public int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, Vector deltaP, double[] du) {
        if (Math.abs(r.getX(0)) > 1) return super.bump(atom, oldState, r, falseTime, deltaP, du);
        Vector v = atom.getVelocity();
        Vector p = atom.getPosition();
        double de = getDeltaU(atom, r, falseTime, oldState);
        double m = atom.getType().getMass();
        if (oldState == 0) {
            // ideal gas trying to become a SQW atom
            if (de < Double.POSITIVE_INFINITY) {
                double v0 = v.getX(0);
                // we (always) have enough kinetic energy to go through the wall
                double newV = Math.sqrt(v0 * v0 - 2 * de / m);
                deltaP.Ea1Tv1(-1, v);
                v.setX(0, newV);
                deltaP.PE(v);
                deltaP.TE(m);
                p.setX(0, p.getX(0) + falseTime * (v0 - newV) + 1e-10);
                du[0] = de;
                updateStatesForSQW(atom);
//                    System.out.println(atom+" IG => SQW "+de+" "+(p.getX(0)+falseTime*v.getX(0)));

                return 1;
            } else {
                // bounce
                du[0] = 0;
                v.setX(0, -v.getX(0));
                deltaP.E(0);
                deltaP.setX(0, 2 * v.getX(0));
                deltaP.TE(m);
                double newP = atom.getPosition().getX(0) - falseTime * v.getX(0) * 2.0;
                p.setX(0, newP - 1e-10);
//                System.out.println(atom+" IG => IG "+(p.getX(0)+falseTime*v.getX(0)));
                return 0;
            }
        } else {
            de = -de;
            // SQW trying to become an ideal gas atom
            double v0 = v.getX(0);
            double ke = 0.5 * m * v0 * v0;
            if (ke > de) {
                // we have enough kinetic energy to go through the wall
                du[0] = de;
                double newV = -Math.sqrt(v0 * v0 - 2 * de / m);
                deltaP.Ea1Tv1(-1, v);
                v.setX(0, newV);
                deltaP.PE(v);
                deltaP.TE(m);
                p.setX(0, p.getX(0) + falseTime * (v0 - newV) - 1e-10);
                clearStatesForIG(atom);
//                System.out.println(atom+" SQW => IG "+de+" "+(p.getX(0)+falseTime*newV));
                return 0;
            } else {
                // bounce
                du[0] = 0;
                v.setX(0, -v.getX(0));
                deltaP.E(0);
                deltaP.setX(0, 2 * v.getX(0));
                deltaP.TE(m);
                double newP = atom.getPosition().getX(0) - falseTime * v.getX(0) * 2.0;
                p.setX(0, newP + 1e-10);
//                System.out.println(atom+" SQW => SQW "+(p.getX(0)+falseTime*v.getX(0)));
                return 1;
            }
        }
    }

    protected void updateStatesForSQW(IAtomKinetic atom) {

        neighborManager.makeNeighborIterator().iterAllNeighbors(atom.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
            @Override
            public void accept(IAtom jAtom, Vector rij) {
                int newState = stateHash.get(jAtom.getLeafIndex());
                if (newState >= 0 && newState < 2) {
                    neighborManager.setPairState(atom.getLeafIndex(), jAtom.getLeafIndex(), newState);
                }
            }
        });
    }

    protected void clearStatesForIG(IAtomKinetic atom) {
        neighborManager.makeNeighborIterator().iterAllNeighbors(atom.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
            @Override
            public void accept(IAtom jAtom, Vector rij) {
                neighborManager.setPairState(atom.getLeafIndex(), jAtom.getLeafIndex(), -1);
            }
        });
    }

    protected double getDeltaU(IAtomKinetic atom, Vector r, double falseTime, int oldState) {

        Boundary boundary = neighborManager.getBox().getBoundary();
        IPotentialPair[] iPotentials = potentials[atom.getType().getIndex()];
        if (oldState == 0) {
            stateHash.clear();
        }
        double de = neighborManager.makeNeighborIterator().iterAndSumAllNeighbors(atom, new NeighborIterator.SuperNbrConsumer() {
            @Override
            public double accept(IAtom atom1, IAtom atom2, Vector rij) {
                Vector dr = Vector.d(r.getD());
                dr.E(atom2.getPosition());
                dr.PEa1Tv1(falseTime, ((IAtomKinetic) atom2).getVelocity());
                if (dr.getX(0) < 0) return 0;
                dr.ME(r);
                boundary.nearestImage(dr);
                P2HardGeneric p2 = (P2HardGeneric) iPotentials[atom2.getType().getIndex()];
                double rCore = p2.getCollisionDiameter(0);
                double rWell = p2.getCollisionDiameter(1);
                double r2 = dr.squared();
                int s = r2 < rCore * rCore ? 0 : (r2 < rWell * rWell ? 1 : -1);
                if (s == 0 && oldState == 1) s = 1;
                double u = p2.getEnergyForState(s);
                if (oldState == 0) {
                    stateHash.put(atom2.getLeafIndex(), s);
                }
                return u;
            }
        });
        return de;
    }
}
