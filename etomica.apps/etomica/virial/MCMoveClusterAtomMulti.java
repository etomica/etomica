package etomica.virial;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;
import etomica.space.IVectorRandom;

/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterAtomMulti extends MCMoveAtom {

    public MCMoveClusterAtomMulti(ISimulation sim, PotentialMaster potentialMaster, int nAtoms) {
        super(sim, potentialMaster);
        selectedAtoms = new IAtomPositioned[nAtoms];
        translationVectors = new IVectorRandom[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            translationVectors[i] = (IVectorRandom)potential.getSpace().makeVector();
        }
        weightMeter = new MeterClusterWeight(potential);
        setStepSize(1.2);
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
            selectedAtoms[i].getPosition().PE(translationVectors[i]);
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
        AtomSet leafList = phase.getSpeciesMaster().getLeafList();
        int total = leafList.getAtomCount();
    	for(int i=1; i<total; i++) {
    		selectedAtoms[i-1] = (IAtomPositioned)leafList.getAtom(i);
    	}
    }

    public void rejectNotify() {
        for(int i=0; i<selectedAtoms.length; i++) {
            selectedAtoms[i].getPosition().ME(translationVectors[i]);
        }
    	((PhaseCluster)phase).rejectNotify();
    }
    
    public void acceptNotify() {
    	((PhaseCluster)phase).acceptNotify();
    }

    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    private final IAtomPositioned[] selectedAtoms;
    private final IVectorRandom[] translationVectors;
}
