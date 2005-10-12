package etomica.virial;

import etomica.integrator.mcmove.MCMoveAtom;
import etomica.phase.Phase;
import etomica.simulation.Simulation;


/**
 *  Overrides MCMoveAtom to prevent index-0 molecule from being displaced
 */

/*
 * Created on Sep 10, 2004
 */

public class MCMoveClusterAtom extends MCMoveAtom implements MCMoveCluster {

    public MCMoveClusterAtom(Simulation sim) {
        super(sim);
        weightMeter = new MeterClusterWeight(sim.potentialMaster);
	}
	
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
    }
    
	public boolean doTrial() {
        PhaseCluster phase = (PhaseCluster)phases[0];
		atom = phase.getSpeciesMaster().atomList.getRandom();
		while(atom.node.getOrdinal()==1) atom = phase.getSpeciesMaster().atomList.getRandom();
		uOld = weightMeter.getDataAsScalar();
        translationVector.setRandomCube();
        translationVector.TE(stepSize);
        atom.coord.position().PE(translationVector);
		phase.trialNotify(atom);
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

    private MeterClusterWeight weightMeter;
}
