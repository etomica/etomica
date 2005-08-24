package etomica.virial;

import etomica.Debug;
import etomica.Phase;
import etomica.atom.AtomSet;
import etomica.potential.Potential0;
import etomica.space.Space;

/**
 * @author David Kofke
 *
 * Pair potential given according to the Mayer bonds in a cluster integral.
 * Does not require that the value of the cluster is non-negative.
 */

/* History
 * 08/20/03 (DAK) small changes to energy method (check for g = 0; abs(g)->g in
 * log argument
 * 08/21/03 (DAK) invoke resetPairs for pairSet in pi method
 * 12/16/03 (DAK) added field to hold indication of whether cluster giving value
 * of potential is positive or negative
 */
public class P0Cluster extends Potential0 {

    private PhaseCluster phaseCluster;
	/**
	 * Constructor for P0Cluster.
	 */
	public P0Cluster(Space space) {
		super(space);
	}
	
    // let's all pretend that the cluster weight is the energy.
	public double energy(AtomSet atoms) {
        if (Debug.ON && atoms.count() > 0) {
            throw new IllegalArgumentException("You actually passed atoms to me.  I'm touched.  Now please stop.");
        }
		return phaseCluster.getSampleCluster().value(phaseCluster.getCPairSet(), phaseCluster.getAPairSet());
	}

    public void setPhase(Phase phase) {
    	phaseCluster = (PhaseCluster)phase;
    }
    
}
