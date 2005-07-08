package etomica.virial.cluster;

import java.util.LinkedList;

import etomica.math.SpecialFunctions;
import etomica.virial.Cluster;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.MayerE;

/**
 * @author kofke
 *
 * Cluster sum that forms the coefficient for an arbitrary coefficient in the
 * activity (z) expansion of the pressure.
 */
public class ClusterBz implements java.io.Serializable {

	/**
	 * Constructor for ClusterBz.
	 * @param weight
	 * @param clusters
	 */
	public static ClusterSum makeClusterBz(int n, MayerE e) {
		ClusterAbstract[] clusters = makeClusters(n,e);
        return new ClusterSum(clusters,weights);
	}

	private static ClusterAbstract[] makeClusters(int j, MayerE e) {
		IntegerSet iSet = new IntegerSet(j);
		int jFact = SpecialFunctions.factorial(j);
		LinkedList list = new LinkedList();
        weights = new double[0];
		do {
			if(!acceptableSet(iSet.set)) continue;
			System.out.println("Setting up "+iSet.toString());
			double currentWeight = (double)coefficient(iSet.set)/(double)jFact;
            double[] weightsTmp = new double[weights.length+1];
            System.arraycopy(weights,0,weightsTmp,0,weights.length);
            weightsTmp[weights.length] = currentWeight;
            weights = weightsTmp;
			Cluster.BondGroup bonds = new Cluster.BondGroup(e, Standard.product(iSet.set));
			Cluster next = new Cluster(j, new Cluster.BondGroup[] {bonds}, true);
			list.add(next);
		} while(iSet.advance());
		int nC = list.size();
		ClusterAbstract[] clusters = new ClusterAbstract[nC];
		clusters = (ClusterAbstract[])list.toArray(clusters);
		return clusters;
	}

    private static double[] weights;
	
	private static boolean acceptableSet(int[] set) {
		int sum = 0;
		final int j = set.length;
		for(int i=1; i<=j; i++) sum += i*set[i-1];
		return sum == j;
	}
	
	private static int coefficient(int[] set) {
		int k = -1;
		final int j = set.length;
		for(int i=1; i<=j; i++) k += set[i-1];
		int sign = ((k % 2)==0) ? +1 : -1; //(-1)^k
		int coeff = sign * SpecialFunctions.factorial(k) * SpecialFunctions.factorial(j);
		for(int i=1; i<=j; i++) {
			int ni = set[i-1];
			if(ni == 0) continue;
			int niFact = SpecialFunctions.factorial(ni);
			int iFactn = power(SpecialFunctions.factorial(i),ni);
			coeff /= (niFact * iFactn);
		}
		return coeff;	
	}
	
	//returns k^n
	private static int power(int k, int n) {
		if(n < 0) throw new IllegalArgumentException();
		int prod = 1;
		for(int i=n; i>0; i--) prod *= k;
		return prod;
	}
		
	
	//class holding a set of integers, with methods to advance 
	//them through all possible values as needed by the cluster
	public static class IntegerSet implements java.io.Serializable {
		public final int n;
		public final int[] set;
		public IntegerSet(int n) {
			this.n = n;
			set = new int[n];
//			for(int i=0; i<n; i++) set[i] = 1;
		}
		
		public String toString() {
			String string = "(";
			for(int i=0; i<n-1; i++) string += set[i]+", ";//System.out.print(set[i]+", ");
			string += set[n-1]+")";//System.out.println(set[n-1]+")");
			return string;
		}
		
		//advances integer set to next value; returns
		//false if already reached end and cannot advance
		public boolean advance() {
			return advance(0);
		}
		//stops advance and moves to next index if (k+1)*set[k] > n
		private boolean advance(int k) {
			if(k >= n) return false;
			set[k]++;
			if((k+1)*set[k]>n) {
				set[k] = 0;
				return advance(k+1);
			} else return true;
		}
		
		//used by BCalculator		
		public boolean advanceFull() {
			return advanceFull(0);
		}
		//stops advance and moves to next index if set[k] > n
		private boolean advanceFull(int k) {
			if(k >= n) return false;
			set[k]++;
			if(set[k]>n) {
				set[k] = 0;
				return advanceFull(k+1);
			} else return true;
		}		
	}//end of IntegerSet
	
	public static void main(String[] args) {
		IntegerSet iSet = new IntegerSet(5);
		do {
			if(acceptableSet(iSet.set)) 
				System.out.println(coefficient(iSet.set)+" "+iSet.toString());
		} while(iSet.advance());
	}
}
