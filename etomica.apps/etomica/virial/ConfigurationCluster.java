package etomica.virial;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomPositionFirstAtom;
import etomica.atom.AtomTypeGroup;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Configuration;
import etomica.config.Conformation;
import etomica.box.Box;
import etomica.space.IVector;
import etomica.space.IVectorRandom;
import etomica.util.IRandom;

/**
 * @author kofke
 *
 * Generates a configuration such that the value of the box's
 * sampling cluster is positive at beta = 1.
 */
public class ConfigurationCluster extends Configuration {

	/**
	 * Constructor for ConfigurationCluster.
	 */
	public ConfigurationCluster(IRandom random) {
		super();
        this.random = random;
	}

	/**
	 * @see etomica.config.Configuration#initializePositions(etomica.AtomIterator)
	 */
    //XXX this can't actually handle multi-atom molecules
	public void initializeCoordinates(Box box) {
		IVectorRandom translationVector = (IVectorRandom)box.getSpace().makeVector();
        IVector dimVector = box.getSpace().makeVector();
        dimVector.E(box.getBoundary().getDimensions());
		IVector center = box.getSpace().makeVector();
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(box);
		iterator.reset();
        for (IAtom a = iterator.nextAtom(); a != null;
             a = iterator.nextAtom()) {
            if (a instanceof IAtomGroup) {
                // initialize coordinates of child atoms
                Conformation config = ((AtomTypeGroup)a.getType()).getConformation();
                config.initializePositions(((IAtomGroup)a).getChildList());
            }
        }
        iterator.reset();

        AtomActionTranslateTo translator = new AtomActionTranslateTo(box.getSpace());
        translator.setDestination(center);
        translator.setAtomPositionDefinition(new AtomPositionFirstAtom());
        translator.actionPerformed(iterator.nextAtom());
        center.E(0.01);
        translator.setDestination(center);
		for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            translator.actionPerformed(atom); //.coord.position().E(center);//put all at center of box
        }
        BoxCluster boxCluster = (BoxCluster)box;
        boxCluster.trialNotify();
     
        if (boxCluster.getSampleCluster() instanceof ClusterWeightAbs) {
    		ClusterAbstract innerCluster = ((ClusterWeightAbs)boxCluster.getSampleCluster()).getWeightCluster();
    		if (innerCluster instanceof ClusterCoupledFlipped) {
    			((ClusterCoupledFlipped)innerCluster).setPhase(boxCluster);
    		}
    }
        
		double value = boxCluster.getSampleCluster().value(boxCluster.getCPairSet(), boxCluster.getAPairSet());
        if (value == 0) {
            System.out.println("initial cluster value bad... trying to fix it.  don't hold your breath.");
        }
		while( value == 0 ) { //if center is not ok, keep trying random positions until ok
		    // if we make it in here we might not make it out!
            iterator.reset();
			iterator.nextAtom();
            for (IAtom a = iterator.nextAtom(); a != null;
                 a = iterator.nextAtom()) {
                translationVector.setRandomCube(random);
                translationVector.TE(dimVector);
                
                translator.setDestination(translationVector);
                translator.actionPerformed(a);
			}
            boxCluster.trialNotify();
			value = boxCluster.getSampleCluster().value(boxCluster.getCPairSet(),boxCluster.getAPairSet());
            System.out.println("value "+value);
            if (value != 0) {
                System.out.println("that wasn't so bad.");
                boxCluster.acceptNotify();
            }
		}
	}

    private static final long serialVersionUID = 3L;
    private final IRandom random;
}
