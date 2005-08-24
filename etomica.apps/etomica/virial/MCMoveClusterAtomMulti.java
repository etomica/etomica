package etomica.virial;

import etomica.Default;
import etomica.atom.Atom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;

/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterAtomMulti extends MCMoveAtom implements MCMoveCluster {

    public MCMoveClusterAtomMulti(PotentialMaster potentialMaster, int nAtoms) {
        super(potentialMaster);
        this.nAtoms = nAtoms;
        selectedAtoms = new Atom[nAtoms];
        translationVectors = new Vector[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            translationVectors[i] = potential.getSpace().makeVector();
        }
        weightMeter = new MeterClusterWeight(potential);
        setStepSize(Default.ATOM_SIZE*1.2);
	}
	
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
		if (selectedAtoms[0] == null) selectAtoms();
        uOld = weightMeter.getDataAsScalar();
        for(int i=0; i<selectedAtoms.length; i++) {
            translationVectors[i].setRandomCube();
            translationVectors[i].TE(stepSize);
            selectedAtoms[i].coord.position().PE(translationVectors[i]);
        }
		((PhaseCluster)phases[0]).trialNotify();
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
    
    public void selectAtoms() {
    	int total=phases[0].getSpeciesMaster().atomList.size();
    	for(int i=total-1; i>total-nAtoms-1; i--) {
    		selectedAtoms[total-1-i] = phases[0].getSpeciesMaster().atomList.get(i);
    	}
    }

    public void rejectNotify() {
        for(int i=0; i<selectedAtoms.length; i++) {
            selectedAtoms[i].coord.position().ME(translationVectors[i]);
        }
    	((PhaseCluster)phases[0]).rejectNotify();
    }
    
    public void acceptNotify() {
    	((PhaseCluster)phases[0]).acceptNotify();
    }

    private MeterClusterWeight weightMeter;
    private final int nAtoms;
    private final Atom[] selectedAtoms;
    private final Vector[] translationVectors;
}
