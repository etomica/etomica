package etomica.virial;

import etomica.Atom;

/**
 * @author kofke
 *
 * Defines the interactions and weight associated with a cluster integral.
 */
public class Cluster {

	/**
	 * Constructs cluster using a single bondgroup, which it wraps in an array
	 * and passes to other constructor.
	 */
	public Cluster(int n, double weight, BondGroup bonds) {
		this(n, weight, new BondGroup[] {bonds});
	}
	/**
	 * Constructor for Cluster.  Each BondGroup in the given array defines a set
	 * of pairs and a single bond function defined between each.
	 */
	public Cluster(int n, double weight, BondGroup[] bonds) {
		super();
		this.n = n;
		this.weight = weight;
		bondArray = new MayerFunction[n][n];
		for(int i=0; i<bonds.length; i++) {
			int[][] idx = bonds[i].pairs;
			if(idx == null) continue;
			for(int j=0; j<idx.length; j++) {
				int i0 = idx[j][0];
				int i1 = idx[j][1];
				if(bondArray[i0][i1]==null) {
					bondArray[i0][i1] = bonds[i].f;
					bondArray[i1][i0] = bonds[i].f;
				} else throw new IllegalArgumentException("Attempting to construct cluster with two bonds defined for a single pair of atoms");
			}
		}
	}
	
	public String toString() {
		String string = "Cluster \n Number of points: " + n + " \n Weight: " + weight +"\n";
		for(int i=0; i<bondArray.length; i++) {
			for(int j=0; j<bondArray.length; j++) {
				string += "("+i+","+j+")"+String.valueOf(bondArray[i][j]) + "\t";
			}
			string += "\n";
		}
		return string;
	}
	public double weight() {return weight;}
	
	/**
	 * Returns the value of the cluster for the given set of atom pairs at the
	 * given value of beta = 1/kT.
	 * @param pairs PairSet defining current configuration.  Does not call reset
	 * for atom pairs.
	 * @param beta reciprocal temperature, 1/kT
	 * @return double value of the cluster
	 */
	public double value(PairSet pairs, double beta) {
		double p = 1.0;
		for(int i=0; i<n-1; i++) {
			for(int j=i+1; j<n; j++) {
				if(bondArray[i][j]==null) continue;
				else p *= bondArray[i][j].f(pairs.getPair(i,j).reset(),beta);
			}
		}
		return p;
	}
	
	/**
	 * Returns the contributions of the given atom for the given set of atom
	 * pairs at the given value of beta = 1/kT.
	 * @param atom Atom for which contribution to cluster is returned by method
	 * @param pairs PairSet defining current configuration.  Does not call reset
	 * for atom pairs.
	 * @param beta reciprocal temperature, 1/kT
	 * @return double value of the cluster
	 */
	public double value(Atom atom, PairSet pairs, double beta) {
		int i = atom.node.index();
		double p = 1.0;
		for(int j=0; j<n; j++) {
			if(bondArray[i][j]==null) continue;
			else p *= bondArray[i][j].f(pairs.getPair(i,j).reset(),beta);
		}
		return p;
	}

	/**
	 * Returns the contributions of the given two atoms for the given set of
	 * atom pairs at the given value of beta = 1/kT.  Returns 1.0 if no bond is
	 * defined between the given pair of atoms for this cluster.
	 * @param atom1 First atom for which contribution to cluster is returned by
	 * method
	 * @param atom2 Second atom for which contribution to cluster is returned by
	 * method
	 * @param pairs PairSet defining current configuration.  Does not call reset
	 * for atom pairs.
	 * @param beta reciprocal temperature, 1/kT
	 * @return double value of the cluster
	 */

	public double value(Atom atom1, Atom atom2, PairSet pairs, double beta) {
		int i = atom1.node.index();
		int j = atom2.node.index();
		return (bondArray[i][j]==null) ? 1.0 : bondArray[i][j].f(pairs.getPair(i,j).reset(),beta);
	}

	private final double weight;
	private final int n;
	private MayerFunction[][] bondArray;
	private PairSet pairs;
	
	/**
	 * Data structure for specifying a bond that is present between one or more
	 * pairs of atoms.
	 * 
	 */
	public static class BondGroup {
		public final MayerFunction f;
		public final int[][] pairs;
		public BondGroup(MayerFunction f, int[][] pairs) {
			this.f = f;
			this.pairs = pairs;
		}
	}
	
}
