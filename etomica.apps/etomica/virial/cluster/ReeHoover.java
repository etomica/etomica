package etomica.virial.cluster;

import etomica.space.CoordinatePair;
import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * Ree-Hoover cluster, in which each pair is joined by either an f-bond or an f-
 * tilde = f+1 bond.
 */
public class ReeHoover extends Cluster {

	/**
	 * Constructor for Ree-Hoover cluster.
	 * @param n number of points in cluster
	 * @param weight weight associated with cluster
	 * @param bonds BondGroup describing the f-bonds in the cluster.  This
	 * constructor will add f-tilde bonds between all pairs not specified with
	 * an f-bond.  The first bondgroup of this cluster will be the given one,
	 * and the second bondgroup will be the f-tilde bonds.
	 */
	public ReeHoover(int n, BondGroup fGroup) {
		super(n, fullBondGroup(n, fGroup));
	}

	/**
	 * Returns a bond group in which every pair not in the given bond group is
	 * given a bond formed from f + 1
	 * @param tilde
	 * @return BondGroup[]
	 */
	private static BondGroup[] fullBondGroup(int n, BondGroup fGroup) {
		int[][] fPairs = fGroup.pairs;
		int nPairs = n*(n-1)/2; //total number of pairs that can be formed from the n points
		int nF = fPairs.length; //number of f bonds that have been specified
		int nTilde = nPairs - nF; //number of f-tilde bonds to be specified
		int[][] fTildePairs = new int[nTilde][];
		int fTildeCount = 0;
		for(int i=0; i<n-1; i++) { //i, j loop over all possible pairs
			for(int j=i+1; j<n; j++) {
				if(notInList(i,j,fPairs)) {//if pair is not in f-bond list, add it to fTilde-bond list
					if(fTildeCount >= nTilde) throw new RuntimeException("Error in Ree-Hoover bond list; perhaps duplicate pair");
					fTildePairs[fTildeCount++] = new int[] {i,j};
				} 
			}
		}
		if(fTildeCount != nTilde) throw new RuntimeException("Error in setup of f-bonds in Ree-Hoover cluster");
		BondGroup fTildeGroup = new BondGroup(new FTilde(fGroup.f), fTildePairs);
		return new BondGroup[] {fGroup, fTildeGroup};
	}//end of fullBondGroup
	
	private static boolean notInList(int i, int j, int[][] pairList) {
		for(int k=0; k<pairList.length; k++) {
			int[] pair = pairList[k];
			if((i==pair[0] && j==pair[1]) || (i==pair[1] && j==pair[0])) return false;//found pair in list
		}
		return true;//didn't find pair in list
	}
	
	private static class FTilde implements MayerFunction {
		private final MayerFunction f;
		private FTilde(MayerFunction f) {
			this.f = f;
		}
		public double f(CoordinatePair cPair, double beta) {
			return f.f(cPair,beta) + 1.0;
		}
		public String toString() {return "f~  ";}
	}//end of FTilde
	
	public static void main(String args[]) {
		Cluster c = new ReeHoover(4, new BondGroup(new etomica.virial.MayerHardSphere(1.0), Standard.ring(4)));
		System.out.println(c.toString());
	}

}
