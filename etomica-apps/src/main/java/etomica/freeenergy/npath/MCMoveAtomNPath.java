package etomica.freeenergy.npath;

import etomica.api.IRandom;
import etomica.atom.AtomSetSinglet;
import etomica.potential.PotentialMaster;
import etomica.space.ISpace;

/**
 * Created by andrew on 4/11/17.
 */
public class MCMoveAtomNPath extends etomica.integrator.mcmove.MCMoveAtom {
    protected final P1ImageHarmonic p1;
    protected final AtomSetSinglet atomSinglet;

    public MCMoveAtomNPath(IRandom random, PotentialMaster potentialMaster, ISpace space, P1ImageHarmonic p1) {
        super(random,potentialMaster,space);
        this.p1 = p1;
        atomSinglet = new AtomSetSinglet();
    }

    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        energyMeter.setTarget(atom);
        atomSinglet.atom = atom;
        uOld = energyMeter.getDataAsScalar()+p1.energy(atomSinglet);
        if(uOld > 1e8 && !fixOverlap) {
            throw new RuntimeException("atom "+atom+" in box "+box+" has an overlap");
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom.getPosition().PE(translationVector);
        return true;
    }//end of doTrial

    public double getB() {
        return super.getB()-p1.energy(atomSinglet);
    }
}
