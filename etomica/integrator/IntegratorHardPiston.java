package etomica.integrator;

import etomica.Atom;
import etomica.Debug;
import etomica.PotentialMaster;
import etomica.potential.P1HardMovingBoundary;
import etomica.potential.PotentialHard;

/**
 * Integrator for DMD with a piston (P1HardMovingBoundary)
 */
public class IntegratorHardPiston extends IntegratorHard {

    /**
     * @param potentialMaster
     * @param potential Potential between piston and every atom in the phase
     */
    public IntegratorHardPiston(PotentialMaster potentialMaster, P1HardMovingBoundary potential) {
        super(potentialMaster);
        pistonPotential = potential;
    }

    public void updateAtom(Atom a) {
        // actually updates the atom
        super.updateAtom(a);
        // check if the atom hit the piston.  if so, then update every atom with the piston
        if (colliderAgent.collisionPotential == pistonPotential) {
            atomIterator.reset();
            listToUpdate.clear();
            while (atomIterator.hasNext()) {
                // look for atoms that wanted to collide with the wall and queue up an uplist recalculation for them.
                Atom[] atom1 = atomIterator.next();
                PotentialHard atom1Potential = ((Agent)atom1[0].ia).collisionPotential;
                if (Debug.ON && Debug.DEBUG_NOW && ((Debug.allAtoms(atom1) && Debug.LEVEL > 1) || (Debug.anyAtom(atom1) && Debug.LEVEL > 2))) {
                    System.out.println(atom1[0]+" thought it would collide with the piston");
                }
                if(atom1Potential == pistonPotential) {
                    if (Debug.ON && Debug.DEBUG_NOW && (Debug.allAtoms(atom1) || Debug.LEVEL > 1)) {
                        System.out.println("Will update "+atom1[0]+" because it wanted to collide with the piston");
                    }
                    listToUpdate.add(atom1[0]);
                }

                
                // recalculate collision time for every atom with the wall
                double collisionTime = pistonPotential.collisionTime(atom1,collisionTimeStep);
                if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atom1)))) {
                    System.out.println("collision down time "+collisionTime+" for atom "+atom1[1]+" with null "+pistonPotential.getClass());
                }
                if(collisionTime < Double.MAX_VALUE) {
                    Agent aia = (Agent)atom1[0].ia;
                    if(collisionTime < aia.collisionTime()) {
                        if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || Debug.anyAtom(atom1))) {
                            System.out.println("setting down time "+collisionTime+" for atom "+atom1[0]+" with null");
                        }
                        if (aia.collisionPotential != null) {
                            aia.eventLinker.remove();
                        }
                        aia.setCollision(collisionTime, atom1[0], pistonPotential);
                        eventList.add(aia.eventLinker);
                    }//end if
                }//end if
            }
            
            processReverseList();
        }
    }
    
    public void advanceAcrossTimeStep(double tStep) {
        super.advanceAcrossTimeStep(tStep);
        pistonPotential.advanceAcrossTimeStep(tStep);
    }
    
    private final P1HardMovingBoundary pistonPotential;
    
}
