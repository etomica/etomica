package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.IOrientation;
import etomica.space.ISpace;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.Orientation3D;

/**
 * Extension of MCMoveAtom that does trial in which several atom orientations are
 * perturbed.  However, orientation of first atom is never altered.  
 */
public class MCMoveClusterAtomRotateMulti extends MCMoveAtom {

    public MCMoveClusterAtomRotateMulti(IRandom random, ISpace _space, int nAtoms) {
        super(null, random, _space, 1.0, Math.PI, false);
        selectedAtoms = new IAtomOriented[nAtoms];
        oldOrientations = new IOrientation[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            oldOrientations[i] = _space.makeOrientation();
        }
        setStepSize(1.2);
	}

    public void setBox(IBox box) {
        super.setBox(box);
        if (selectedAtoms[0] == null) selectAtoms();
        for(int i=0; i<selectedAtoms.length; i++) {
            if (selectedAtoms[i].getOrientation() instanceof Orientation3D) {
                oldOrientations[i] = new Orientation3D(space);
            }
        }
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
		if (selectedAtoms[0] == null) selectAtoms();
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        for(int i=0; i<selectedAtoms.length; i++) {
            oldOrientations[i].E(selectedAtoms[i].getOrientation());
            selectedAtoms[i].getOrientation().randomRotation(random, stepSize);
        }
        if (random.nextInt(100) == 0) {
            for(int i=0; i<selectedAtoms.length; i++) {
                ((IVectorMutable)selectedAtoms[i].getOrientation().getDirection()).normalize();
            }
        }
		((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
		return true;
	}
	
    public double getA() {
        return uNew/uOld;
    }

    public double getB() {
    	return 0.0;
    }
    
    public void selectAtoms() {
        IAtomList leafList = box.getLeafList();
        int total = leafList.getAtomCount();
    	for(int i=0; i<total; i++) {
    		selectedAtoms[i] = (IAtomOriented)leafList.getAtom(i);
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

    private final IAtomOriented[] selectedAtoms;
    private final IOrientation[] oldOrientations;
}
