package etomica.virial;

import etomica.api.IAtomSet;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.ISimulation;

import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.PotentialMaster;


/**
 *  Overrides MCMoveAtom to prevent index-0 molecule from being displaced
 */
public class MCMoveClusterAtom extends MCMoveAtom {

    public MCMoveClusterAtom(ISimulation sim, PotentialMaster potentialMaster) {
        super(sim, potentialMaster);
        weightMeter = new MeterClusterWeight(potentialMaster);
	}
	
    public void setBox(IBox p) {
        super.setBox(p);
        weightMeter.setBox(p);
    }
    
	public boolean doTrial() {
        IAtomSet leafList = box.getLeafList();
		atom = leafList.getAtom(random.nextInt(1+leafList.getAtomCount()-1));
		uOld = weightMeter.getDataAsScalar();
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        ((IAtomPositioned)atom).getPosition().PE(translationVector);
		((BoxCluster)box).trialNotify();
		uNew = Double.NaN;
		return true;
	}
	
    public double getB() {
    	return 0;
    }
    
    public double getA() {
        uNew = weightMeter.getDataAsScalar();
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }

    public void rejectNotify() {
    	super.rejectNotify();
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	super.acceptNotify();
    	((BoxCluster)box).acceptNotify();
    }

    private static final long serialVersionUID = 1L;
    private MeterClusterWeight weightMeter;
}
