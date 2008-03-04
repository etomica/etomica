package etomica.virial;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IRandom;

import etomica.atom.IAtomOriented;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.PotentialMaster;
import etomica.space.IOrientation;

/**
 * Extension of MCMoveAtom that does trial in which several atom orientations are
 * perturbed.  However, orientation of first atom is never altered.  
 */
public class MCMoveClusterAtomRotateMulti extends MCMoveAtom {

    public MCMoveClusterAtomRotateMulti(IRandom random, PotentialMaster potentialMaster, int nAtoms) {
        super(potentialMaster, random, 1.0, Math.PI, false);
        selectedAtoms = new IAtomOriented[nAtoms];
        oldOrientations = new IOrientation[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            oldOrientations[i] = potential.getSpace().makeOrientation();
        }
        weightMeter = new MeterClusterWeight(potential);
        setStepSize(1.2);
	}
	
    public void setBox(IBox p) {
        super.setBox(p);
        weightMeter.setBox(p);
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
		if (selectedAtoms[0] == null) selectAtoms();
        uOld = weightMeter.getDataAsScalar();
        for(int i=0; i<selectedAtoms.length; i++) {
            oldOrientations[i].E(selectedAtoms[i].getOrientation());
            selectedAtoms[i].getOrientation().randomRotation(random, stepSize);
        }
		((BoxCluster)box).trialNotify();
        uNew = weightMeter.getDataAsScalar();
		return true;
	}
	
    public double getA() {
        return uNew/uOld;
    }

    public double getB() {
    	return 0.0;
    }
    
    public void selectAtoms() {
        IAtomSet leafList = box.getLeafList();
        int total = leafList.getAtomCount();
    	for(int i=1; i<total; i++) {
    		selectedAtoms[i-1] = (IAtomOriented)leafList.getAtom(i);
    	}
    }

    public void rejectNotify() {
        for(int i=0; i<selectedAtoms.length; i++) {
            selectedAtoms[i].getOrientation().E(oldOrientations[i]);
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    private final IAtomOriented[] selectedAtoms;
    private final IOrientation[] oldOrientations;
}
