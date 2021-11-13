/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorHard;
import etomica.potential.IPotentialHard;
import etomica.potential.IPotentialHardField;
import etomica.potential.compute.NeighborManagerHard;
import etomica.space.Vector;
import etomica.species.SpeciesManager;
import etomica.util.random.IRandom;

/**
 * Integrator for DMD with a piston (P1HardMovingBoundary)
 */
public class IntegratorHardPiston extends IntegratorHard {

    public IntegratorHardPiston(IPotentialHard[][] pairPotentials, IPotentialHardField[] fieldPotentials,
                                NeighborManagerHard neighborManager, IRandom random, double timeStep, double temperature, Box box,
                                etomica.modules.pistoncylinder.P1HardMovingBoundary pistonPotential, SpeciesManager sm) {
        super(pairPotentials, fieldPotentials, neighborManager, random, timeStep, temperature, box, sm, null);
        setTimeStep(1.0);
        this.pistonPotential = pistonPotential;
    }

    public void resetPiston() {
        if (space.D() == 3) {
            pistonPotential.setWallPosition(box.getBoundary().getBoxSize().getX(1) * 0.5);
        } else {
            pistonPotential.setWallPosition(-box.getBoundary().getBoxSize().getX(1) * 0.5);
        }
    }

    public void doStepInternal() {
        if (pistonUpdateRequested) {
            pistonUpdateRequested = false;
            updatePiston(0);
        }
        super.doStepInternal();
    }

    @Override
    protected void updateAtom(IAtomKinetic a, double falseTime) {
        boolean doUpdatePiston = collisionPartners[a.getLeafIndex()] < 0 && !pistonPotential.isStationary();
        // actually updates the atom
        super.updateAtom(a, falseTime);
        // check if the atom hit the piston.  if so, then update every atom with the piston
        if (doUpdatePiston) {
            updatePiston(falseTime);
        }
    }

    public void reset() {
        super.reset();
        updatePiston(0);
    }

    /**
     * recalculate collision times for all atoms with the wall/piston
     */
    public void updatePiston(double falseTime) {
        AtomArrayList downColliders = new AtomArrayList(0);
        // look for atoms that wanted to collide with the wall and queue up an uplist recalculation for them.
        IAtomList leafList = box.getLeafList();
        Vector r = box.getSpace().makeVector();
        int nAtoms = leafList.size();
        for (int iLeaf = 0; iLeaf < nAtoms; iLeaf++) {
            IAtomKinetic atom1 = (IAtomKinetic) leafList.get(iLeaf);

            IPotentialHard atom1Potential = collisionPotentials[iLeaf];
            if (atom1Potential == null) {
                downColliders.add(atom1);
            }

            // recalculate collision time for every atom with the wall
            r.E(atom1.getPosition());
            r.PEa1Tv1(falseTime, atom1.getVelocity());
            double collisionTime = pistonPotential.collisionTime(atom1, r, atom1.getVelocity(), 0, falseTime);
            if (collisionTime < 0) {
                System.out.println("computed piston update " + iLeaf + " at t=" + collisionTime + " at r=" + r);
                pistonPotential.collisionTime(atom1, r, atom1.getVelocity(), 0, falseTime);
                throw new RuntimeException("negative!");
            }

            double tcol = collisionTime + falseTime + tBase;
            if (tcol < Math.min(collisionTimes[iLeaf], tMax)) {
                if (collisionTimes[iLeaf] < tMax) {
                    // need to remove j from event bins
                    removeCollision(iLeaf);
                }
                collisionTimes[iLeaf] = tcol;
                collisionPartners[iLeaf] = -1;
                collisionPotentials[iLeaf] = null;
                collisionOldState[iLeaf] = 0;
                r.PEa1Tv1(collisionTime, atom1.getVelocity());
                collisionVector[iLeaf].E(r);

                int b = (int) (collisionTimes[iLeaf] / tMax * eventBinsFirstAtom.length);

                int oldFirst = eventBinsFirstAtom[b];
                eventBinsNextAtom[iLeaf] = oldFirst;
                eventBinsPrevAtom[iLeaf] = -1 - b;
                if (oldFirst >= 0) eventBinsPrevAtom[oldFirst] = iLeaf;
                eventBinsFirstAtom[b] = iLeaf;
            }
        }

        for (IAtom iAtom : downColliders) {
            collisionTimeUp((IAtomKinetic) iAtom, falseTime);
        }
    }

    public void advanceAcrossTimeStep(double tStep) {
        super.advanceAcrossTimeStep(tStep);
        pistonPotential.advanceAcrossTimeStep(tStep);
    }

    public synchronized void pistonUpdateRequested() {
        pistonUpdateRequested = true;
    }

    private final P1HardMovingBoundary pistonPotential;
    private boolean pistonUpdateRequested = false;
}
