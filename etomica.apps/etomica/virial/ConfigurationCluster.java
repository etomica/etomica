package etomica.virial;

import etomica.Atom;
import etomica.Configuration;
import etomica.Space;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListCompound;
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
        iterator = new AtomIteratorListCompound();
	}

	/**
	 * @see etomica.Configuration#initializePositions(etomica.AtomIterator)
	 */
    //XXX this can't actually handle multi-atom molecules
	public void initializePositions(AtomList[] lists) {
		Vector translationVector = phase.space().makeVector();
        Vector dimVector = Space.makeVector(dimensions);
		Vector center = phase.space().makeVector();
		iterator.setLists(lists);
		iterator.reset();
        if (!iterator.hasNext()) return;
		while(iterator.hasNext()) iterator.nextAtom().coord.position().E(center);//put all at center of box
        phase.trialNotify();
		double value = phase.getSampleCluster().value(phase.getCPairSet(), 1.0);
        if (value == 0) {
            System.out.println("initial cluster value bad... trying to fix it.  don't hold your breath.");
        }
		while( value == 0 ) { //if center is not ok, keep trying random positions until ok
		    // if we make it in here we probably won't make it out!
            iterator.reset();
			iterator.nextAtom();
			while(iterator.hasNext()) {
                translationVector.setRandomCube();
                translationVector.TE(dimVector);
                Atom a = iterator.nextAtom();
                a.coord.position().PE(translationVector);
//                System.out.println("moved "+a+" to "+a.coord.position());
			}
			value = phase.getSampleCluster().value(phase.getCPairSet(),1.0);
            if (value != 0) {
                System.out.println("that wasn't so bad.");
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
    private final AtomIteratorListCompound iterator;
}
