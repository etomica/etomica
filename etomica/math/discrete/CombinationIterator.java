package etomica.math.discrete;

import etomica.math.SpecialFunctions;
/**
 * @author kofke
 *
 * Iterator that returns different combinations of n integers taken k at a time.
 * For example, if constructed with n = 4 and k = 2, will return, in successive
 * calls to next(), the arrays {0,1},{0,2},{0,3},{1,2},{1,3},{2.3}.  Entries in
 * each array will be in numerical order.
 */

/* History
 * 10/30/03 (DAK) new
 */
 
public class CombinationIterator {

	private int[] index;
	private int n, k;
		
	/**
	 * Constructor for CombinationIterator.
	 */
	public CombinationIterator() {
		super();
	}

	public void reset(int n, int k) {
		if(n < 0 || k > n) throw new IllegalArgumentException("CombinationIterator must have 0 <= k <= n");
		this.n = n;
		this.k = k;
		index = new int[k];
		for(int i=0; i<k; i++) index[i] = i;
		if(k != 0) index[k-1]--;
	}
	
	public boolean hasNext() {
		if(index == null) return false;
		if(k == 0) return true;
		for(int j=0; j<k; j++) if(index[j] < n-k+j) return true;
		return false;
	}
	
	public int[] next() {
		if(k == 0) {
			index = null;
			return new int[] {};
		}
		int j = k-1;
		while(index[j] == n-k+j) j--;
		index[j]++;
		for(int i=j+1; i<k; i++) index[i] = index[i-1]+1;
		return index;
	}
	
	public static void main(String[] args) {
//		int n = Integer.parseInt(args[0]);
//		int k = Integer.parseInt(args[1]);
		int n = 10;
		int k = 10;
		CombinationIterator p = new CombinationIterator();
		p.reset(n,k);
		int count = 1;
		while(p.hasNext()) {
			int[] a = p.next();
			System.out.print(count++ +". {");
			for(int i=0; i<a.length-1; i++) System.out.print(a[i]+", ");
			if(a.length > 0) System.out.print(a[a.length-1]);
			System.out.println("}");
		}
		int s = SpecialFunctions.factorial(n)/SpecialFunctions.factorial(k)/SpecialFunctions.factorial(n-k);
		System.out.println("Expected total:" + s);
	}

}
