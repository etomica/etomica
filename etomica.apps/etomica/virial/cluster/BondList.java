package etomica.virial.cluster;

import etomica.math.discrete.CombinationIterator;
import etomica.virial.*;
import etomica.virial.cluster.graphics.*;

/**
 * Data structure that holds list of atom pairs that are joined by bonds.
 * Also has several utility methods to charactize a cluster formed from these
 * bonds.
 */
public class BondList {
	
	private int[][] pairs;

	/**
	 * Constructor for BondList.
	 */
	public BondList(int[][] pairs) {
		super();
		this.pairs = pairs;
	}

	public int bondCount() {
		return pairs.length;
	}
	
	public boolean isPermutationOf(BondList test, int[] rootPoints) {
		return false;
	}

	/**
	 * Returns the set of all n-point clusters in the h-bond expansion of the
	 * bridge function.
	 * @param n
	 * @return Cluster[]
	 */
	private static Cluster[] methodName(int n) {
		int[][] temp = Standard.full(n);//starting diagram -- all others formed by deleting bonds from it
		int[][] allBonds = new int[temp.length-1][];//but first start by deleting 0-1 bond (between root points)
		for(int i=0, k=0; i<temp.length; i++) if(temp[i][0]+temp[i][1]!=1) allBonds[k++] = temp[i];//eliminates {0,1} as these are the only ones that sum to 1

		DeletionIterator iterator = new DeletionIterator();
		iterator.reset(new BondList(allBonds),allBonds.length-(n-1));
		
		java.util.LinkedList list = new java.util.LinkedList();
		while(iterator.hasNext()) {//loop over all graphs obtained by deleted from full graph
			int[][] bonds = iterator.next();
			if(!isConnected(n, bonds)) continue;//eliminate unconnected diagrams
			java.util.Iterator listIterator = list.iterator();
			boolean unique = true;
//			while(listIterator.hasNext()) {//loop over stored permutation sets
//				PermutationSet perm = (PermutationSet)listIterator.next();
//				if(perm.contains(bonds)) {unique = false; break;} 
//			}
//			if(unique) list.add(new PermutationSet(bonds));
		}
		//linked list has all unique clusters at this point
		//now check each for meeting restrictions of h-bond elementary diagrams
		
		//bridge 
		java.util.Iterator listIterator = list.iterator();
		java.util.LinkedList acceptable = new java.util.LinkedList();
		while(listIterator.hasNext()) {
//			int[][] bonds = ((PermutationSet)listIterator.next()).bonds();
////			if(hasNode(bonds)) continue;
//			if(hasArticulationPoint(n, bonds)) continue;
//			if(hasArticulationPairs(n, bonds)) continue;
//			acceptable.add(bonds);
		}
		
		listIterator = acceptable.iterator();
		Cluster[] clusters = new Cluster[acceptable.size()];
		int i = 0;
		return null;
//		while(listIterator.hasNext()) clusters[i++] = new Cluster(n, 1.0, new Cluster.BondGroup(null, (int[][])listIterator.next()));
	}
/**
 * Iterator that generates bond lists beginning from an initial list and
 * deleting one or more bonds from it.
 */
	private static class DeletionIterator {
		
		private CombinationIterator comboIterator = new CombinationIterator();
		private BondList next;
		private int bondCount, originalBondCount;
		private int nDelete;
		private int maxDelete;
		private int[][] bonds;
		private int[][] allBonds;
		private boolean hasNext;
		
		public int[][] next() {
			int[] deletedBonds = comboIterator.next();
			int k = 0;//index for bonds being kept
			int j = 0;//pointer to deleted bonds
			if(bonds.length != bondCount) bonds = new int[bondCount][];
			for(int i=0; i<bondCount; i++) {//copy nondeleted pairs to array for next
				while((j<nDelete) && k == deletedBonds[j]) {k++; j++;}//loop over bonds being deleted until finding one to keep
				bonds[i] = allBonds[k++];//copy bond that isn't being deleted
			}
			if(!comboIterator.hasNext()) {
				if(nDelete == maxDelete) hasNext = false;
				else comboIterator.reset(originalBondCount, ++nDelete);
				bondCount = originalBondCount - nDelete;
			}
			return bonds;
		}
		
		public void reset(BondList seed, int maxDelete) {
			allBonds = seed.pairs;
			originalBondCount = seed.bondCount();
			this.maxDelete = maxDelete;
			nDelete = 0;
			comboIterator.reset(originalBondCount, nDelete);
			bondCount = originalBondCount - nDelete;
			bonds = new int[bondCount][];	
			hasNext = true;		
		}
		
		public boolean hasNext() {
			return hasNext;
		}
		
	}//end of DeletionIterator
	
	/**
	 * Examines array of bonds to see if they form a connected graph.
	 * @param n number of points subject to bonds
	 * @param bonds array of bonds between points
	 * @return boolean true if bonds described a connected graph
	 */
	public static boolean isConnected(int n, int[][] bonds) {
		//make an connectivity array
		boolean[][] connected = new boolean[n][n];
		for(int i=0; i<bonds.length; i++) {
			int j0 = bonds[i][0];
			int j1 = bonds[i][1];
			connected[j0][j1] = connected[j1][j0] = true;
		}
		//loop to see if we can get a connection between point 0 and all other points
		int oldCount = -1;
		int count = 0;
		while(oldCount != count && (count < (n-1))) {//end loop if updating of 0's connections results in no change
			oldCount = count;
			count = 0;
			for(int i=1; i<n; i++) {
				if(connected[0][i]) {//0 is connected to i; update array to show that 0 is connected to everything i is connected to
					for(int j=0; j<n; j++) connected[0][j] |= connected[i][j];//this is an "or =" operator
					count++;
				} 
			}
		}
		return (count == (n-1));//true if 0 is connected to all others by some path		
	}
	
	public static class PermutationSet {
		final int[][][] signatures;
		final int[][] prototype;
		final int[][] nbrCount; //nbrCount[i][0] is no of root point i is bonded to, [i][1] is no of field points
		final int nBonds, n;
//		java.util.LinkedList permutations = new java.util.LinkedList();
		PermutationSet(int n, int[][] bonds) {
			this.n = n;
			nBonds = bonds.length;
			prototype = new int[nBonds][];
			for(int i=0; i<nBonds; i++) {
				prototype[i][0] = bonds[i][0];
				prototype[i][1] = bonds[i][1];
			}
			int[][] nbrCount = new int[n][2];
			for(int i=0; i<nBonds; i++) {//loop over bonds
				int p1 = bonds[i][0];//label of each point in bond
				int p2 = bonds[i][1];
				//increment count of points for bonds to root (=0,1) points or field points
				if(p2 == 0 || p2 == 1) nbrCount[p1][0]++; else nbrCount[p1][1]++;
				if(p1 == 0 || p1 == 1) nbrCount[p2][0]++; else nbrCount[p2][1]++;
			}
			int[][][] signatures = new int[n][][];
			for(int j=0; j<n; j++) {
				
			}
			//make all permutations
		}
		
		public boolean contains(int[][] test) {
			java.util.Iterator iterator = permutations.iterator();
			while(iterator.hasNext()) {
				int[][] next = (int[][])iterator.next();
				if(equal(test,next)) return true;
			}
			return false;
		}
	}
	
//	public static void main(String[] args) {
//		DeletionIterator iterator = new DeletionIterator();
//		BondList bonds = new BondList(Standard.full(4));
//		iterator.reset(bonds, 6);
//		while(iterator.hasNext()) {
//			int[][] a = iterator.next();
//			System.out.print(isConnected(4,a)+", {");
//			for(int i=0; i<a.length-1; i++)  System.out.print("{"+a[i][0]+","+a[i][1]+"}, ");
//			if(a.length > 0) System.out.print("{"+a[a.length-1][0]+","+a[a.length-1][1]+"}");
//			System.out.println("}");
//		}
//	}
		public static void main(String[] args) {
			int nPoints = 4;
			int maxDelete = 6;
			java.util.LinkedList testList = new java.util.LinkedList();

			BondList bonds = new BondList(Standard.full(4));
			DeletionIterator iterator = new DeletionIterator();
			iterator.reset(bonds, maxDelete);
			while(iterator.hasNext()) {
				int[][] a = iterator.next();

//				System.out.print(isConnected(4,a)+ ", {");
//				for(int i=0; i<a.length-1; i++)  System.out.print("{"+a[i][0]+","+a[i][1]+"}, ");
//				if(a.length > 0) System.out.print("{"+a[a.length-1][0]+","+a[a.length-1][1]+"}");
//				System.out.println("}");

				int[][] b = new int[a.length][2];
				for(int i=0;i<a.length;i++){
					b[i][0]=a[i][0];				
					b[i][1]=a[i][1];
				}
				testList.add(b);
			}			
	
			ShowCluster showCluster = new ShowCluster(nPoints, testList);		
				showCluster.addWindowListener(new java.awt.event.WindowAdapter(){
				public void windowClosing( java.awt.event.WindowEvent e){
					System.exit(0);
				}});
		}// end of main

}
