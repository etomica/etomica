package etomica.virial;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.ISimulation;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterAtomMulti extends MCMoveAtom {

    public MCMoveClusterAtomMulti(ISimulation sim, IPotentialMaster potentialMaster,
    		                      ISpace _space) {
        super(sim, potentialMaster, _space);
        weightMeter = new MeterClusterWeight(potential);
        setStepSize(1.2);
	}
	
    public void setBox(IBox p) {
        super.setBox(p);
        weightMeter.setBox(p);
        translationVectors = new IVectorRandom[box.getMoleculeList().getAtomCount()-1];
        for (int i=0; i<box.getMoleculeList().getAtomCount()-1; i++) {
            translationVectors[i] = (IVectorRandom)space.makeVector();
        }
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
        uOld = weightMeter.getDataAsScalar();
        IAtomList leafAtoms = box.getLeafList();
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
        IAtomList leafAtoms = box.getLeafList();
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
