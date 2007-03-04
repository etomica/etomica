package etomica.virial;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.phase.Phase;
import etomica.simulation.Simulation;


/**
 *  Overrides MCMoveAtom to prevent index-0 molecule from being displaced
 */

/*
 * Created on Sep 10, 2004
 */

public class MCMoveClusterAtom extends MCMoveAtom {

    public MCMoveClusterAtom(Simulation sim) {
        super(sim);
        weightMeter = new MeterClusterWeight(sim.getPotentialMaster());
	}
	
    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
    }
    
	public boolean doTrial() {
        AtomArrayList leafList = phase.getSpeciesMaster().leafList;
		atom = leafList.get(random.nextInt(1+leafList.size()-1));
		uOld = weightMeter.getDataAsScalar();
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        ((AtomLeaf)atom).getCoord().getPosition().PE(translationVector);
		((PhaseCluster)phase).trialNotify();
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
    	((PhaseCluster)phase).rejectNotify();
    }
    
    public void acceptNotify() {
    	super.acceptNotify();
    	((PhaseCluster)phase).acceptNotify();
    }

    private static final long serialVersionUID = 1L;
    private MeterClusterWeight weightMeter;
}
