package etomica.virial;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.IVectorRandom;

/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterAtomMulti extends MCMoveAtom {

    public MCMoveClusterAtomMulti(Simulation sim, int nAtoms) {
        super(sim);
        selectedAtoms = new AtomLeaf[nAtoms];
        translationVectors = new IVectorRandom[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            translationVectors[i] = (IVectorRandom)potential.getSpace().makeVector();
        }
        weightMeter = new MeterClusterWeight(potential);
        setStepSize(sim.getDefaults().atomSize*1.2);
	}
	
    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
		if (selectedAtoms[0] == null) selectAtoms();
        uOld = weightMeter.getDataAsScalar();
        for(int i=0; i<selectedAtoms.length; i++) {
            translationVectors[i].setRandomCube(random);
            translationVectors[i].TE(stepSize);
            selectedAtoms[i].getCoord().getPosition().PE(translationVectors[i]);
        }
		((PhaseCluster)phase).trialNotify();
		uNew = Double.NaN;
		return true;
	}
	
    public double getA() {
        uNew = weightMeter.getDataAsScalar();
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }

    public double getB() {
    	return 0.0;
    }
    
    public void selectAtoms() {
        AtomArrayList leafList = phase.getSpeciesMaster().getLeafList();
        int total = leafList.size();
    	for(int i=1; i<total; i++) {
    		selectedAtoms[i-1] = (AtomLeaf)leafList.get(i);
    	}
    }

    public void rejectNotify() {
        for(int i=0; i<selectedAtoms.length; i++) {
            selectedAtoms[i].getCoord().getPosition().ME(translationVectors[i]);
        }
    	((PhaseCluster)phase).rejectNotify();
    }
    
    public void acceptNotify() {
    	((PhaseCluster)phase).acceptNotify();
    }

    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    private final AtomLeaf[] selectedAtoms;
    private final IVectorRandom[] translationVectors;
}
