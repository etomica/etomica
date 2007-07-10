package etomica.virial;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;


/**
 *  Overrides MCMoveAtom to prevent index-0 molecule from being displaced
 */
public class MCMoveClusterAtom extends MCMoveAtom {

    public MCMoveClusterAtom(ISimulation sim, PotentialMaster potentialMaster) {
        super(sim, potentialMaster);
        weightMeter = new MeterClusterWeight(potentialMaster);
	}
	
    public void setBox(Box p) {
        super.setBox(p);
        weightMeter.setBox(p);
    }
    
	public boolean doTrial() {
        AtomSet leafList = box.getLeafList();
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
