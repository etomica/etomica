package etomica.virial;

/**
 * @author kofke
 *
 * Defines the interactions and weight associated with a cluster integral.
 */
public class Cluster {

	/**
	 * Constructor for Cluster.
	 */
	public Cluster(double weight, int[][] pairs) {
		super();
		this.weight = weight;
		this.pairs = pairs;
		nPairs = pairs.length;
	}
	
	public double weight() {return weight;}
	
	public double value(double[][] f) {
		double p = 1.0;
		for(int i=0; i<nPairs; i++) {
			int[] pair = pairs[i];
			int i0 = pair[0];
			int i1 = pair[1]-i0-1;
			p *= f[i0][i1];
		}
		return p;
	}	

	private final double weight;
	private final int[][] pairs;
	private final int nPairs;
}
