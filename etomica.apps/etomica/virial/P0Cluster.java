package etomica.virial;

import etomica.api.IAtomSet;
import etomica.api.IBox;

import etomica.potential.Potential0;
import etomica.space.Space;

/**
 * @author David Kofke
 *
 * Pair potential given according to the Mayer bonds in a cluster integral.
 * Does not require that the value of the cluster is non-negative.
 */
public class P0Cluster extends Potential0 {

    private static final long serialVersionUID = 1L;
    private BoxCluster boxCluster;
	/**
	 * Constructor for P0Cluster.
	 */
	public P0Cluster(Space space) {
		super(space);
	}
	
    // let's all pretend that the cluster weight is the energy.
	public double energy(IAtomSet atoms) {
        return 0;
	}

    public double weight() {
        return boxCluster.getSampleCluster().value(boxCluster);
    }

    public void setBox(IBox box) {
    	boxCluster = (BoxCluster)box;
    }
}
