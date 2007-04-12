package etomica.virial;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomGroup;
import etomica.atom.AtomPositionFirstAtom;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Configuration;
import etomica.config.Conformation;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.IVectorRandom;
import etomica.util.IRandom;

/**
 * @author kofke
 *
 * Generates a configuration such that the value of the phase's
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
	public void initializeCoordinates(Phase phase) {
		IVectorRandom translationVector = (IVectorRandom)phase.getSpace().makeVector();
        IVector dimVector = phase.getSpace().makeVector();
        dimVector.E(phase.getBoundary().getDimensions());
		IVector center = phase.getSpace().makeVector();
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(phase);
		iterator.reset();
        while (iterator.hasNext()) {
            IAtom a = iterator.nextAtom();
            if (!a.isLeaf()) {
                // initialize coordinates of child atoms
                Conformation config = a.getType().creator().getConformation();
                config.initializePositions(((AtomGroup)a).getChildList());
            }
        }
        iterator.reset();

        AtomActionTranslateTo translator = new AtomActionTranslateTo(phase.getSpace());
        translator.setDestination(center);
        translator.setAtomPositionDefinition(new AtomPositionFirstAtom());
        translator.actionPerformed(iterator.nextAtom());
        center.E(0.01);
        translator.setDestination(center);
        if (!iterator.hasNext()) return;
		while(iterator.hasNext()) { 
            translator.actionPerformed(iterator.nextAtom()); //.coord.position().E(center);//put all at center of box
        }
        PhaseCluster phaseCluster = (PhaseCluster)phase;
        phaseCluster.trialNotify();
		double value = phaseCluster.getSampleCluster().value(phaseCluster.getCPairSet(), phaseCluster.getAPairSet());
        if (value == 0) {
            System.out.println("initial cluster value bad... trying to fix it.  don't hold your breath.");
        }
		while( value == 0 ) { //if center is not ok, keep trying random positions until ok
		    // if we make it in here we might not make it out!
            iterator.reset();
			iterator.nextAtom();
			while(iterator.hasNext()) {
                translationVector.setRandomCube(random);
                translationVector.TE(dimVector);
                IAtom a = iterator.nextAtom();
                
                translator.setDestination(translationVector);
                translator.actionPerformed(a);
			}
            phaseCluster.trialNotify();
			value = phaseCluster.getSampleCluster().value(phaseCluster.getCPairSet(),phaseCluster.getAPairSet());
            System.out.println("value "+value);
            if (value != 0) {
                System.out.println("that wasn't so bad.");
                phaseCluster.acceptNotify();
            }
		}
	}

    private static final long serialVersionUID = 3L;
    private final IRandom random;
}
