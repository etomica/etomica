package etomica.integrator.mcmove;

import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.atom.IAtomOriented;
import etomica.space.IOrientation;
import etomica.space.ISpace;

/**
 * Performs a rotation of an atom (not a molecule) that has an orientation coordinate.
 */
public class MCMoveRotate extends MCMoveAtom {
    
    private static final long serialVersionUID = 2L;
    private final IOrientation oldOrientation;

    private transient IOrientation iOrientation;

    public MCMoveRotate(IPotentialMaster potentialMaster, IRandom random,
    		            ISpace _space) {
        super(potentialMaster, random, _space, Math.PI/2, Math.PI, false);
        oldOrientation = _space.makeOrientation();
    }
    
    public boolean doTrial() {
        if(box.getMoleculeList().getMoleculeCount()==0) {return false;}
        atom = atomSource.getAtom();

        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        iOrientation = ((IAtomOriented)atom).getOrientation(); 
        oldOrientation.E(iOrientation);  //save old orientation
        iOrientation.randomRotation(random, stepSize);
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public void rejectNotify() {
        iOrientation.E(oldOrientation);
    }
}
