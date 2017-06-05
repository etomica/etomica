package etomica.integrator.mcmove;

import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * MC move that operates only on atoms within a given region of x values.
 * If the selected atom is moved outside the region, it is wrapped back inside.
 * 
 * @author Andrew Schultz
 */
public class MCMoveAtomInRegion extends MCMoveAtom {

    protected int nAttempts;
    protected double xMin, xMax;
    protected Vector oldPosition;
    
    public MCMoveAtomInRegion(IRandom random, PotentialMaster potentialMaster,
            Space _space) {
        super(random, potentialMaster, _space);
        oldPosition = _space.makeVector();
    }

    /**
     * Defines the region atoms must be selected from.  The atom is chosen
     * randomly.  If the atom is not in the region, repeated attempts are
     * made until one is found or maxAttempts is reached.
     * 
     * If {@code newXMax > newXMin}, then {@code newXMax > x > newXMin}.
     * If {@code newXMax < newXMin}, then {@code newXMax < x < newXMin} (the region spans the
     * boundary, excluding the middle of the box)
     */
    public void setXRange(double newXMin, double newXMax, int maxAttempts) {
        xMin = newXMin;
        xMax = newXMax;
        nAttempts = maxAttempts;
    }
    
    public boolean doTrial() {
        boolean success = false;
        for (int i=0; i<nAttempts; i++) {
            atom = atomSource.getAtom();
            if (atom == null) return false;
            double x = atom.getPosition().getX(0);
            if (xMin < xMax) {
                if (x > xMin && x < xMax) {
                    success = true;
                    break;
                }
            }
            else {
                if (x < xMax || x > xMin) {
                    success = true;
                    break;
                }
            }
        }
        if (!success) return false;
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        if(uOld > 1e8 && !fixOverlap) {
            throw new RuntimeException("atom "+atom+" in box "+box+" has an overlap");
        }
        oldPosition.E(atom.getPosition());
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom.getPosition().PE(translationVector);
        Vector dx = box.getBoundary().centralImage(atom.getPosition());
        atom.getPosition().PE(dx);
        double newX = atom.getPosition().getX(0);
        if (xMin < xMax) {
            if (newX < xMin) {
                // too negative
                newX = xMax - (xMin-newX);
            }
            else if (newX > xMax) {
                // too positive
                newX = xMin + (newX-xMax);
            }
        }
        else if (newX < xMin && newX > xMax) {
            // in the middle (excluded)
            double newX1 = xMax - (xMin-newX);
            double newX2 = xMin + (newX-xMax);
            newX = (-newX1 > newX2) ? newX2 : newX1;
        }
        atom.getPosition().setX(0,newX);
        translationVector.Ev1Mv2(atom.getPosition(), oldPosition);
        
        return true;
    }

}
