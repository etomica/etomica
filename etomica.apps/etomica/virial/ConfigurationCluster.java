package etomica.virial;

import etomica.AtomIterator;
import etomica.Configuration;
import etomica.Space;
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
	}

	/**
	 * @see etomica.Configuration#initializePositions(etomica.AtomIterator)
	 */
	public void initializePositions(AtomIterator[] iterator) {
		Vector translationVector = phase.space().makeVector();
        Vector dimVector = Space.makeVector(dimensions);
		Vector center = phase.space().makeVector();
		AtomIterator iter = iterator[0];
		iter.reset();
		while(iter.hasNext()) iter.nextAtom().coord.position().E(center);//put all at center of box
		double value = phase.getSampleCluster().value(phase.getCPairSet(), 1.0);
        if (value == 0) {
            System.out.println("initial cluster value bad... trying to fix it.  don't hold your breath.");
        }
		while( value == 0 ) { //if center is not ok, keep trying random positions until ok
		    // if we make it in here we probably won't make it out!
            iter.reset();
			iter.nextAtom();
			while(iter.hasNext()) {
                translationVector.setRandomCube();
                translationVector.TE(dimVector);
                iter.nextAtom().coord.position().PE(translationVector);
			}
			value = phase.getSampleCluster().value(phase.getCPairSet(),1.0);
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
}
