package etomica.virial;

import etomica.Atom;
import etomica.math.discrete.PermutationIterator;

/**
 * @author kofke
 *
 * Defines the interactions associated with a cluster integral.
 */

/* History
 * 08/19/03 (DAK) added capability for using permutations
 * 08/21/03 (DAK) removed calls to reset() for pairs when used in value method.
 * This requires that the class calling these method must do the reset first.
 * Change is made for efficiency, as in some cases the several cluster values
 * are computed for the same configuration.  Calling resetPairs once before
 * computing them all thus saves on recalculation of the r2 separations.
 */
 
public class ClusterBonds implements java.io.Serializable {

	/**
	 * Constructor for Cluster, with default to not using permutations (this can
	 * be changed after construction if desired).
	 */
	public ClusterBonds(int n, int[][][] bonds) {
		this(n, bonds, false);
	}
	
	/**
	 * Constructor for Cluster.  Each BondGroup in the given array defines a set
	 * of pairs and a single bond function defined between each.  
	 */	
	public ClusterBonds(int n, int[][][] bonds, boolean usePermutations) {
		super();
		this.nPoints = n;
		bondIndexArray = new int[n][n];
		for(int i=0; i<n; i++) for(int j=0; j<n; j++) bondIndexArray[i][j] = -1;//index value that indicates no bond between pair
		nBondTypes = bonds.length;
		for(int i=0; i<nBondTypes; i++) {
			int[][] idx = bonds[i];
			if(idx == null) continue;
			for(int j=0; j<idx.length; j++) {
				int i0 = idx[j][0];
				int i1 = idx[j][1];
				bondIndexArray[i0][i1] = i;
				bondIndexArray[i1][i0] = i;
			}
		}
		if(usePermutations) setUsePermutations(usePermutations);
	}
	
	public String toString() {
		String string = "Cluster \n Number of points: " + nPoints + "\n";
		if(!usePermutations) {
			for(int i=0; i<nPoints; i++) {
				for(int j=0; j<nPoints; j++) {
                    int protoIndex = bondIndexArray[i][j];
                    if (protoIndex > -1) {
                        string += "("+i+","+j+") "+protoIndex + "\t";
                    }
				}
				string += "\n";
			}
		} else {
			string += "Using permutations \n";
			for(int s=0; s<nPermutations; s++) {
				int[] p = permutations[s];// e.g. {0,2,1,4} 
				for(int i=0; i<nPoints; i++) {
					int ip = p[i];//index of atom in position i for permuted labeling (e.g., i=1, ip=2}
					for(int j=0; j<nPoints; j++) {
						int jp = p[j];
						int protoIndex = bondIndexArray[i][j];//index of bond between i and j in prototype diagram
						if (protoIndex > -1) {
						    string += "("+ip+","+jp+") "+ protoIndex  + "\t";
                        }
					}
					string += "\n";
				}//end for i
				string += "\n";
			}//end for s
		}//end else
		return string;
	}
	
	
	public int pointCount() {return nPoints;}
	
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
	public double value(double[][][] fValues) {
		if(usePermutations) return valueUsingPermutations(fValues);
		double p = 1.0;
		for(int i=0; i<nPoints-1; i++) {
			for(int j=i+1; j<nPoints; j++) {
                int protoIndex = bondIndexArray[i][j];
			    if(protoIndex == -1) continue;
			    p *= fValues[i][j][protoIndex];
			    if(p == 0.0) return 0.0;
			}
		}
		return p;
	}
	
	/**
	 * Returns value obtained by averaging over all unique permutations of the
	 * point indices.
	 * @param pairs
	 * @param beta
	 * @return double
	 */
	private double valueUsingPermutations(double[][][] fValues) {
		double sum = 0;
		
		//loop over permutations
		for(int s=0; s<nPermutations; s++) {
			int[] p = permutations[s];// e.g. {0,2,1,4} 
			double prod = 1.0;
			pairLoop: for(int i=0; i<nPoints-1; i++) {
				double[][] fi = fValues[p[i]];//p[i] = index of atom in position i for permuted labeling (e.g., i=1, ip=2}
				int[] bi = bondIndexArray[i];
				for(int j=i+1; j<nPoints; j++) {
					int protoIndex = bi[j];//index of bond between i and j in prototype diagram
//                    System.out.println(i+" "+j+" protoIndex="+protoIndex);
					if(protoIndex < 0) continue;//no bond there
//                    System.out.println("v="+fi[p[j]][protoIndex]);
					prod *= fi[p[j]][protoIndex];//value of bond for permuted diagram
					if(prod == 0.0) break pairLoop;
				}
			}//end pairLoop
			sum += prod;
/*            if (prod != 0.0 && (bondGroup[0].f instanceof MayerHardSphere || bondGroup[0].f instanceof MayerEHardSphere)) {
                System.out.println("whoo-hoo!!!!!!!!!!!!!  we got one!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                throw new RuntimeException();
            }*/
		}
		sum *= rPermutations;//divide by nPermutations
		return sum;
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
	public double value(Atom atom, double[][][] fPairs) {
		int i = atom.node.index();
		double p = 1.0;
		for(int j=0; j<nPoints; j++) {
            int protoIndex = bondIndexArray[i][j];
			if(protoIndex==-1) continue;
			p *= fPairs[protoIndex][i][j];
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

	public double value(Atom atom1, Atom atom2, double[][][] fPairs) {
		int i = atom1.node.index();
		int j = atom2.node.index();
        int protoIndex = bondIndexArray[i][j];
		return (protoIndex==-1) ? 1.0 : fPairs[protoIndex][i][j];
	}

	private final int nPoints; //number of points (molecules) in cluster
 	protected final int[][] bondIndexArray;//array giving bondGroup index of each bond in bondArray
					                       //bondArray[i][j] = bondGroup[bondIndexArray[i][j]]
					                       //bondIndexArray[i][j] < 0 indicates bondArray[i][j] == null
	private final int nBondTypes; //length of bondGroup (usually 1, maybe 2)
	private boolean usePermutations = false;//flag indicating if value is via prototype cluster only, or by average of all its permutations
	private int[][] permutations;//array of permutations of atom indexes giving cluster different from prototype
	private int nPermutations; //permutations.length
	private double rPermutations; //reciprocal of nPermutations
	
	private void makePermutations() {
		java.util.LinkedList pList = new java.util.LinkedList();//permutations giving unique arrays
		java.util.LinkedList aList = new java.util.LinkedList();//unique permuted arrays
		
		PermutationIterator pIter = new PermutationIterator(nPoints);//iterator over permutations
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
		rPermutations = 1/(double)nPermutations;
	}//end of makePermutations

	//returns true if all elements in two arrays are equal term by term
	private boolean intArrayEqual(int[][] a0, int[][] a1) {
		for(int i=0; i<nPoints; i++) for(int j=0; j<nPoints; j++) if(a0[i][j] != a1[i][j]) return false;
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
//		Cluster cluster = new etomica.virial.cluster.Ring(5, new MayerHardSphere(1.0));
//		Cluster cluster = new etomica.virial.cluster.ReeHoover(4, new BondGroup(new MayerHardSphere(1.0),etomica.virial.cluster.Standard.chain(4)));
//		cluster.setUsePermutations(true);
//		System.out.println(cluster.toString());
	}
}
