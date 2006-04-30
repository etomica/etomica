package etomica.virial;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPositionFirstAtom;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorArrayListCompound;
import etomica.config.Configuration;
import etomica.config.Conformation;
import etomica.space.Space;
import etomica.space.Vector;

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
	public ConfigurationCluster(Space space) {
		super(space);
        iterator = new AtomIteratorArrayListCompound();
	}

	/**
	 * @see etomica.config.Configuration#initializePositions(etomica.AtomIterator)
	 */
    //XXX this can't actually handle multi-atom molecules
	public void initializePositions(AtomArrayList[] lists) {
		Vector translationVector = phase.space().makeVector();
        Vector dimVector = Space.makeVector(dimensions);
		Vector center = phase.space().makeVector();
		iterator.setLists(lists);
		iterator.reset();
        while (iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            if (!a.node.isLeaf()) {
                // initialize coordinates of child atoms
                Conformation config = a.type.creator().getConformation();
                config.initializePositions(((AtomTreeNodeGroup) a.node).childList);
            }
        }
        iterator.reset();

        AtomActionTranslateTo translator = new AtomActionTranslateTo(space);
        translator.setDestination(center);
        translator.setAtomPositionDefinition(new AtomPositionFirstAtom());
        translator.actionPerformed(iterator.nextAtom());
        center.E(0.01);
        translator.setDestination(center);
        if (!iterator.hasNext()) return;
		while(iterator.hasNext()) { 
            translator.actionPerformed(iterator.nextAtom()); //.coord.position().E(center);//put all at center of box
        }
        phase.trialNotify();
		double value = phase.getSampleCluster().value(phase.getCPairSet(), phase.getAPairSet());
        if (value == 0) {
            System.out.println("initial cluster value bad... trying to fix it.  don't hold your breath.");
        }
		while( value == 0 ) { //if center is not ok, keep trying random positions until ok
		    // if we make it in here we might not make it out!
            iterator.reset();
			iterator.nextAtom();
			while(iterator.hasNext()) {
                translationVector.setRandomCube();
                translationVector.TE(dimVector);
                Atom a = iterator.nextAtom();
                
                translator.setDestination(translationVector);
                translator.actionPerformed(a);
			}
            phase.trialNotify();
			value = phase.getSampleCluster().value(phase.getCPairSet(),phase.getAPairSet());
            System.out.println("value "+value);
            if (value != 0) {
                System.out.println("that wasn't so bad.");
                phase.acceptNotify();
            }
		}
	}

	/**
	 * Returns the phase.
	 * @return Phase
	 */
	public PhaseCluster getPhase() {
		return phase;
	}

    /**
	 * Sets the phase.
	 * @param phase The phase to set
	 */
	public void setPhase(PhaseCluster phase) {
		this.phase = phase;
	}

    private PhaseCluster phase;
    private final AtomIteratorArrayListCompound iterator;
}
