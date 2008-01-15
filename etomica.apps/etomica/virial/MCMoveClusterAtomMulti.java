package etomica.virial;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.box.Box;
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

    public MCMoveClusterAtomMulti(ISimulation sim, PotentialMaster potentialMaster) {
        super(sim, potentialMaster);
        weightMeter = new MeterClusterWeight(potential);
        setStepSize(1.2);
	}
	
    public void setBox(Box p) {
        super.setBox(p);
        weightMeter.setBox(p);
        translationVectors = new IVectorRandom[box.getMoleculeList().getAtomCount()-1];
        for (int i=0; i<box.getMoleculeList().getAtomCount()-1; i++) {
            translationVectors[i] = (IVectorRandom)potential.getSpace().makeVector();
        }
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
        uOld = weightMeter.getDataAsScalar();
        AtomSet leafAtoms = box.getLeafList();
        for(int i=1; i<leafAtoms.getAtomCount(); i++) {
            translationVectors[i-1].setRandomCube(random);
            translationVectors[i-1].TE(stepSize);
            ((IAtomPositioned)leafAtoms.getAtom(i)).getPosition().PE(translationVectors[i-1]);
        }
		((BoxCluster)box).trialNotify();
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
    
    public void rejectNotify() {
        AtomSet leafAtoms = box.getLeafList();
        for(int i=1; i<leafAtoms.getAtomCount(); i++) {
            ((IAtomPositioned)leafAtoms.getAtom(i)).getPosition().ME(translationVectors[i-1]);
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    private IVectorRandom[] translationVectors;
}
