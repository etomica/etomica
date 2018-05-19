/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardMovingBoundary;
import etomica.potential.PotentialHard;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.util.Debug;

/**
 * Integrator for DMD with a piston (P1HardMovingBoundary)
 */
public class IntegratorHardPiston extends IntegratorHard {

    /**
     * @param potentialMaster
     * @param potentialMaster Potential between piston and every atom in the box
     */
    public IntegratorHardPiston(Simulation sim,
                                PotentialMaster potentialMaster,
                                P1HardMovingBoundary potentialWrapper, Box box) {
        super(sim, potentialMaster, box);
        setTimeStep(1.0);
        pistonPotential = potentialWrapper;
        atomSetSinglet = new AtomSetSinglet();
    }

    public void resetPiston() {
        if (space.D() == 3) {
            pistonPotential.setWallPosition(box.getBoundary().getBoxSize().getX(1)*0.5);
        }
        else {
            pistonPotential.setWallPosition(-box.getBoundary().getBoxSize().getX(1)*0.5);
        }
    }

    protected void doStepInternal() {
        if (pistonUpdateRequested) {
            pistonUpdateRequested = false;
            updatePiston();
        }
        super.doStepInternal();
    }
    
    public void updateAtom(IAtom a) {
        boolean isPistonPotential = agentManager.getAgent(a).collisionPotential == pistonPotential;
        // actually updates the atom
        super.updateAtom(a);
        // check if the atom hit the piston.  if so, then update every atom with the piston
        if (isPistonPotential) {
            updatePiston();
        }
    }
    
    public void reset() {
        super.reset();
        updatePiston();
    }
    
    /**
     * recalculate collision times for all atoms with the wall/piston
     */
    public void updatePiston() {
        listToUpdate.clear();
        // look for atoms that wanted to collide with the wall and queue up an uplist recalculation for them.
        IAtomList leafList = box.getLeafList();
        int nAtoms = leafList.size();
        for (int iLeaf=0; iLeaf<nAtoms; iLeaf++) {
            IAtom atom1 = leafList.get(iLeaf);
            atomSetSinglet.atom = atom1;
            PotentialHard atom1Potential = agentManager.getAgent(atom1).collisionPotential;
            if (Debug.ON && Debug.DEBUG_NOW && ((Debug.allAtoms(atomSetSinglet) && Debug.LEVEL > 1) || (Debug.anyAtom(atomSetSinglet) && Debug.LEVEL > 2))) {
                System.out.println(atom1+" thought it would collide with the piston");
            }
            if(atom1Potential == pistonPotential) {
                if (Debug.ON && Debug.DEBUG_NOW && (Debug.allAtoms(new AtomSetSinglet(atom1)) || Debug.LEVEL > 1)) {
                    System.out.println("Will update "+atom1+" because it wanted to collide with the piston");
                }
                listToUpdate.add(atom1);
            }

            
            // recalculate collision time for every atom with the wall
            double collisionTime = pistonPotential.collisionTime(atomSetSinglet,collisionTimeStep);
            if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atomSetSinglet)))) {
                System.out.println("collision down time "+collisionTime+" for atom "+atom1+" with null "+pistonPotential.getClass());
            }
            if(collisionTime < Double.POSITIVE_INFINITY) {
                Agent aia = agentManager.getAgent(atom1);
                if(collisionTime < aia.collisionTime()) {
                    if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || Debug.anyAtom(atomSetSinglet))) {
                        System.out.println("setting down time "+collisionTime+" for atom "+atom1+" with null");
                    }
                    if (aia.collisionPotential != null) {
                        aia.eventLinker.remove();
                    }
                    aia.setCollision(collisionTime, null, pistonPotential);
                    eventList.add(aia.eventLinker);
                }//end if
            }//end if
        }
        
        processReverseList();
    }
    
    public void advanceAcrossTimeStep(double tStep) {
        super.advanceAcrossTimeStep(tStep);
        pistonPotential.advanceAcrossTimeStep(tStep);
    }
    
    public synchronized void pistonUpdateRequested() {
        pistonUpdateRequested = true;
    }
    
    private static final long serialVersionUID = 1L;
    private final P1HardMovingBoundary pistonPotential;
    private boolean pistonUpdateRequested = false;
    protected final AtomSetSinglet atomSetSinglet;
}
