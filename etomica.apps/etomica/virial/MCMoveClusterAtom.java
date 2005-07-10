package etomica.virial;

import etomica.Phase;
import etomica.PotentialMaster;
import etomica.integrator.mcmove.MCMoveAtom;


/**
 *  Overrides MCMoveAtom to prevent index-0 molecule from being displaced
 */

/*
 * Created on Sep 10, 2004
 */

public class MCMoveClusterAtom extends MCMoveAtom implements MCMoveCluster {
	private MeterClusterWeight weightMeter;
	public MCMoveClusterAtom(PotentialMaster potentialMaster) {
		super(potentialMaster);
        weightMeter = new MeterClusterWeight(potentialMaster);
	}
	
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
    }
    
	public boolean doTrial() {
        Phase phase = phases[0];
		atom = phase.getSpeciesMaster().atomList.getRandom();
		while(atom.node.index()==0) atom = phase.getSpeciesMaster().atomList.getRandom();
		// this slows things down due to caching
//		weightMeter.setTarget(atom);
		uOld = weightMeter.getDataAsScalar();
        translationVector.setRandomCube();
        translationVector.TE(stepSize);
        atom.coord.position().PE(translationVector);
		((PhaseCluster)phase).trialNotify();
		uNew = Double.NaN;
		return true;
	}
	
    public double lnProbabilityRatio() {
    	return Math.log(probabilityRatio());
    }
    
    public double probabilityRatio() {
        uNew = weightMeter.getDataAsScalar();
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }

    public double trialRatio() {
    	return 1.0;
    }

    public void rejectNotify() {
    	super.rejectNotify();
    	((PhaseCluster)phases[0]).rejectNotify();
    }
    
    public void acceptNotify() {
    	super.acceptNotify();
    	((PhaseCluster)phases[0]).acceptNotify();
    }
}
