package etomica.virial;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.Configuration;
import etomica.Default;
import etomica.Simulation;
import etomica.Space;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class ConfigurationCluster extends Configuration {

	/**
	 * Constructor for ConfigurationCluster.
	 * @param space
	 */
	public ConfigurationCluster(Space space) {
		super(space);
	}

	/**
	 * Constructor for ConfigurationCluster.
	 * @param sim
	 */
	public ConfigurationCluster(Simulation sim) {
		super(sim);
	}

	/**
	 * @see etomica.Configuration#initializePositions(etomica.AtomIterator)
	 */
	public void initializePositions(AtomIterator[] iterator) {
		Space.Vector dimensions = phase.space.makeVector();
		dimensions.E(Default.BOX_SIZE);
		Space.Vector center = phase.space.makeVector();
		center.Ea1Tv1(0.5, dimensions);
		AtomIterator iter = iterator[0];
		iter.reset();
		while(iter.hasNext()) iter.next().coord.translateTo(center);//put all at center of box
		double value = cluster.value(phase.getPairSet(), 1.0);
		while( value == 0 || (signPositive != (value>0.0))) { //if center is not ok, keep trying random positions until ok
			iter.reset();
			iter.next().coord.translateTo(center);
			while(iter.hasNext()) {
				Atom a = iter.next();
				a.coord.translateToRandom(phase);	
			}
			value = cluster.value(phase.getPairSet(),1.0);
		}//end while
	}

	boolean signPositive;
	Cluster cluster;
	PhaseCluster phase;
	/**
	 * Returns the cluster.
	 * @return Cluster
	 */
	public Cluster getCluster() {
		return cluster;
	}

	/**
	 * Returns the phase.
	 * @return Phase
	 */
	public PhaseCluster getPhase() {
		return phase;
	}

	/**
	 * Returns the signPositive.
	 * @return boolean
	 */
	public boolean isSignPositive() {
		return signPositive;
	}

	/**
	 * Sets the cluster.
	 * @param cluster The cluster to set
	 */
	public void setCluster(Cluster cluster) {
		this.cluster = cluster;
	}

	/**
	 * Sets the phase.
	 * @param phase The phase to set
	 */
	public void setPhase(PhaseCluster phase) {
		this.phase = phase;
	}

	/**
	 * Sets the signPositive.
	 * @param signPositive The signPositive to set
	 */
	public void setSignPositive(boolean signPositive) {
		this.signPositive = signPositive;
	}

}
