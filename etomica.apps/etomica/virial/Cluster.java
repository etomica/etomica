package etomica.virial;

import etomica.Atom;
import etomica.math.discrete.PermutationIterator;

/**
 * @author kofke
 *
 * Defines the interactions and weight associated with a cluster integral.
 */

/* History
 * 08/19/03 (DAK) added capability for using permutations
 * 08/21/03 (DAK) removed calls to reset() for pairs when used in value method.
 * This requires that the class calling these method must do the reset first.
 * Change is made for efficiency, as in some cases the several cluster values
 * are computed for the same configuration.  Calling resetPairs once before
 * computing them all thus saves on recalculation of the r2 separations.
 */
 
public class Cluster implements ClusterAbstract {

	/**
	 * Constructs cluster using a single bondgroup, which it wraps in an array
	 * and passes to other constructor.
	 */
	public Cluster(int n, double weight, BondGroup bonds) {
		this(n, weight, new BondGroup[] {bonds});
	}
	/**
	 * Constructor for Cluster, with default to not using permutations (this can
	 * be changed after construction if desired).
	 */
	public Cluster(int n, double weight, BondGroup[] bonds) {
		this(n, weight, bonds, false);
	}
	
	/**
	 * Constructor for Cluster.  Each BondGroup in the given array defines a set
	 * of pairs and a single bond function defined between each.  
	 */	
	public Cluster(int n, double weight, BondGroup[] bonds, boolean usePermutations) {
		super();
		this.n = n;
		this.weight = weight;
		bondCount = 0;
		for(int i=0; i<bonds.length; i++) bondCount += bonds[i].pairs.length;
		bondArray = new MayerFunction[n][n];
		bondIndexArray = new int[n][n];
		for(int i=0; i<n; i++) for(int j=0; j<n; j++) bondIndexArray[i][j] = -1;//index value that indicates no bond between pair
		nBondTypes = bonds.length;
		bondGroup = new BondGroup[nBondTypes];
		for(int i=0; i<nBondTypes; i++) {
			bondGroup[i] = new BondGroup(bonds[i]);
			bondGroup[i].index = i;
			int[][] idx = bonds[i].pairs;
			if(idx == null) continue;
			for(int j=0; j<idx.length; j++) {
				int i0 = idx[j][0];
				int i1 = idx[j][1];
				if(bondArray[i0][i1]==null) {
					bondArray[i0][i1] = bonds[i].f;
					bondArray[i1][i0] = bonds[i].f;
					bondIndexArray[i0][i1] = i;
					bondIndexArray[i1][i0] = i;
				} else throw new IllegalArgumentException("Attempting to construct cluster with two bonds defined for a single pair of atoms");
			}
		}
		if(usePermutations) setUsePermutations(usePermutations);
	}
	
	/**
	 * Makes this cluster using the given Mayer function joining the pairs given
	 * by the first bondgroup of the given cluster.
	 * @param f Mayer function for this cluster
	 * @param cluster Cluster defining bonds to which the Mayer function is
	 * applied.
	 */
	public Cluster(MayerFunction f, Cluster cluster) {
		this(cluster.n, cluster.weight, new BondGroup(f, cluster.bondGroup[0].pairs));
	}
	
	public String toString() {
		String string = "Cluster \n Number of points: " + n + " \n Weight: " + weight +"\n";
		if(!usePermutations) {
			for(int i=0; i<n; i++) {
				for(int j=0; j<n; j++) {
					string += "("+i+","+j+")"+String.valueOf(bondArray[i][j]) + "\t";
				}
				string += "\n";
			}
		} else {
			string += "Using permutations \n";
			for(int s=0; s<nPermutations; s++) {
				int[] p = permutations[s];// e.g. {0,2,1,4} 
				for(int i=0; i<n; i++) {
					int ip = p[i];//index of atom in position i for permuted labeling (e.g., i=1, ip=2}
					for(int j=0; j<n; j++) {
						int jp = p[j];
						int protoIndex = bondIndexArray[i][j];//index of bond between i and j in prototype diagram
						MayerFunction bond = (protoIndex < 0) ? null : bondGroup[protoIndex].f;
						string += "("+ip+","+jp+")"+String.valueOf(bond) + "\t";
					}
					string += "\n";
				}//end for i
				string += "\n";
			}//end for s
		}//end else
		return string;
	}
	
	
	public double weight() {return weight;}
	
	public int pointCount() {return n;}
	
	public int bondCount() {return bondCount;}
	
	public BondGroup[] bondGroup() {return bondGroup;}
	
	public boolean hasOddBondCount() {return (bondCount % 2) != 0;}
		
	/**
	 * Returns the value of the cluster for the given set of atom pairs at the
	 * given value of beta = 1/kT. It is expected that the pairs in the given
	 * PairSet will have all be reset (i.e., their r2 values already computed).
	 * @param pairs PairSet defining current configuration.  Does not call reset
	 * for atom pairs, so this must be done before calling the method
	 * (accomplished by invoking pairs.resetPairs() method).
	 * @param beta reciprocal temperature, 1/kT
	 * @return double value of the cluster
	 */
	public double value(PairSet pairs, double beta) {
		if(usePermutations) return valueUsingPermutations(pairs, beta);
		else {
			double p = 1.0;
			for(int i=0; i<n-1; i++) {
				for(int j=i+1; j<n; j++) {
					if(bondArray[i][j]==null) continue;
					else p *= bondArray[i][j].f(pairs.getPair(i,j),beta);
					
					if(p == 0.0) return 0.0;
				}
			}
			return p;
		}
	}
	
	/**
	 * Returns value obtained by averaging over all unique permutations of the
	 * point indices.
	 * @param pairs
	 * @param beta
	 * @return double
	 */
	private double valueUsingPermutations(PairSet pairs, double beta) {
		double sum = 0;
		
		//evaluate and store values of all bonds between all pairs
		for(int i=0; i<n-1; i++) {
			for(int j=i+1; j<n; j++) {
				for(int k=0; k<nBondTypes; k++) {
					f[i][j][k] = bondGroup[k].f.f(pairs.getPair(i,j),beta);
					f[j][i][k] = f[i][j][k];
				}
			}
		}
		
		//loop over permutations
		for(int s=0; s<nPermutations; s++) {
			int[] p = permutations[s];// e.g. {0,2,1,4} 
			double prod = 1.0;
			pairLoop: for(int i=0; i<n-1; i++) {
				double[][] fj = f[p[i]];//p[i] = index of atom in position i for permuted labeling (e.g., i=1, ip=2}
				int[] bi = bondIndexArray[i];
				for(int j=i+1; j<n; j++) {
					int protoIndex = bi[j];//index of bond between i and j in prototype diagram
					if(protoIndex < 0) continue;//no bond there
					else prod *= fj[p[j]][protoIndex];//value of bond for permuted diagram
					if(prod == 0.0) break pairLoop;
				}
			}//end pairLoop
			sum += prod;
		}
		return sum*rPermutations;//divide by nPermutations
	}//end valueUsingPermutations
	
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
			else p *= bondArray[i][j].f(pairs.getPair(i,j),beta);
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
		return (bondArray[i][j]==null) ? 1.0 : bondArray[i][j].f(pairs.getPair(i,j),beta);
	}

	private final double weight; //weight assigned to cluster
	private final int n; //number of points (molecules) in cluster
	protected int bondCount;//number of bonds in cluster (not to be confused with number of types of bonds, nBondTypes)
 	protected final MayerFunction[][] bondArray;//array of bonds between pairs in prototype (unpermuted) cluster 
 	protected final int[][] bondIndexArray;//array giving bondGroup index of each bond in bondArray
					                       //bondArray[i][j] = bondGroup[bondIndexArray[i][j]]
					                       //bondIndexArray[i][j] < 0 indicates bondArray[i][j] == null
	private BondGroup[] bondGroup;//array with info about all bonds, their index, and pairs they apply to
	private final int nBondTypes; //length of bondGroup (usually 1, maybe 2)
	private double[][][] f;//stores value of every bond between every pair; used by valueUsingPermutations method
	private boolean usePermutations = false;//flag indicating if value is via prototype cluster only, or by average of all its permutations
	private int[][] permutations;//array of permutations of atom indexes giving cluster different from prototype
	private int nPermutations; //permutations.length
	private double rPermutations; //reciprocal of nPermutations
	
	/**
	 * Data structure for specifying a bond that is present between one or more
	 * pairs of atoms.
	 */
	public static class BondGroup {
		public final MayerFunction f;
		public final int[][] pairs;
		public int index;
		//copy constructor
		public BondGroup(BondGroup original) {
			this(original.f, original.pairs);
		}
		public BondGroup(MayerFunction f, int[][] pairs) {
			this.f = f;
			this.pairs = pairs;
		}
		public static final BondGroup NULL = new BondGroup(new MayerHardSphere(0.0), new int[][] {});
	}
	/**
	 * Checks for equality of this cluster with another.  Returns true if they
	 * have the same number of points, the same weight, and if all bonds between
	 * pairs are equal.  Not a test of topological equality (i.e., equality
	 * without considering labels or indexes of points); each point is
	 * considered labeled and clusters must have corresponding bonds between
	 * pairs to return true.
	 */
	public boolean equals(Object obj) {
		if(super.equals(obj)) return true;//return true if this is the same instance as the test cluster
		if(!(obj instanceof Cluster)) return false;
		Cluster test = (Cluster)obj;
		if(this.n != test.n) return false;
		if(this.weight != test.weight) return false;
		for(int i=0; i<n-1; i++) {
			for(int j=i+1; j<n; j++) {
				MayerFunction thisBond = this.bondArray[i][j];
				MayerFunction testBond = test.bondArray[i][j];
				if((thisBond==null) ? (testBond==null) : thisBond.equals(testBond)) continue;//both bonds are null, or they're equal; keep going
				else return false; //mismatch found
			}
		}
		return true;
		//could use java.util.Arrays.equals to test bondArray equality
	}
	
	private void makePermutations() {
		java.util.LinkedList pList = new java.util.LinkedList();//permutations giving unique arrays
		java.util.LinkedList aList = new java.util.LinkedList();//unique permuted arrays
		
		PermutationIterator pIter = new PermutationIterator(n);//iterator over permutations
		pLoop: while(pIter.hasNext()) {//loop over permutations
			int[] p = pIter.next();
			int[][] pIntBondArray = PermutationIterator.matrixPermutation(p, bondIndexArray);//permute array according to permutations
			java.util.Iterator iter = aList.iterator();
			while(iter.hasNext()) {//check new array against previous unique ones
				if(intArrayEqual((int[][])iter.next(), pIntBondArray)) continue pLoop;
			}
			pList.add(p);//permutation is different from prototype; keep it
			aList.add(pIntBondArray); //save array for checking uniqueness of subsequent arrays
		}
		nPermutations = pList.size();
		permutations = new int[nPermutations][];
		pList.toArray(permutations);//fill permutations array with the elements of pList
		rPermutations = 1.0/(double)nPermutations;
		f = new double[n][n][nBondTypes];
	}//end of makePermutations

	//returns true if all elements in two arrays are equal term by term
	private boolean intArrayEqual(int[][] a0, int[][] a1) {
		for(int i=0; i<n; i++) for(int j=0; j<n; j++) if(a0[i][j] != a1[i][j]) return false;
		return true;
	}
	/**
	 * Returns the usePermutations flag.
	 * @return boolean
	 */
	public boolean isUsePermutations() {
		return usePermutations;
	}

	/**
	 * Sets the usePermutations flag.
	 * @param usePermutations The usePermutations to set
	 */
	public void setUsePermutations(boolean usePermutations) {
		this.usePermutations = usePermutations;
		if(usePermutations) makePermutations();
		else nPermutations = 1;
	}
	
	/**
	 * Number of permutations of (and including) the basic cluster used when
	 * calculating its value for a configuration
	 * @return int the number of permutations
	 */
	public int nPermutations() {return nPermutations;}

	public static void main(String[] args) {
		
//		Cluster cluster = new etomica.virial.cluster.D6(new MayerHardSphere(1.0));
//		Cluster cluster = new etomica.virial.cluster.Ring(5, 1.0, new MayerHardSphere(1.0));
		Cluster cluster = new etomica.virial.cluster.ReeHoover(4, 1.0, new BondGroup(new MayerHardSphere(1.0),etomica.virial.cluster.Standard.chain(4)));
		cluster.setUsePermutations(true);
		System.out.println(cluster.toString());
	}
}
