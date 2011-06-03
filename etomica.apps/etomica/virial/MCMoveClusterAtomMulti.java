package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IRandom;
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

    public MCMoveClusterAtomMulti(IRandom random, ISpace _space) {
        super(random, null, _space);
        setStepSize(1.2);
	}
	
    public void setBox(IBox p) {
        super.setBox(p);
        if (translationVectors == null) {
            translationVectors = new IVectorRandom[box.getLeafList().getAtomCount()-1];
            for (int i=0; i<translationVectors.length; i++) {
                translationVectors[i] = (IVectorRandom)space.makeVector();
            }
        }
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IAtomList leafAtoms = box.getLeafList();
        for(int i=1; i<leafAtoms.getAtomCount(); i++) {
            translationVectors[i-1].setRandomCube(random);
            translationVectors[i-1].TE(stepSize);
            leafAtoms.getAtom(i).getPosition().PE(translationVectors[i-1]);
        }
		((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
		return true;
	}
	
    public double getA() {
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }

    public double getB() {
    	return 0.0;
    }
    
    public void rejectNotify() {
        IAtomList leafAtoms = box.getLeafList();
        for(int i=1; i<leafAtoms.getAtomCount(); i++) {
            leafAtoms.getAtom(i).getPosition().ME(translationVectors[i-1]);
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    private static final long serialVersionUID = 1L;
    protected IVectorRandom[] translationVectors;
}
