package etomica.math.discrete;

/**
 * @author kofke
 *
 * Iterator that returns different permutations of a sequence of integers when
 * called in successive iterations.  For example, if constructed with n = 3,
 * will return, in successive calls to next(), the arrays {0,1,2}, {0,2,1},
 * {2,0,1}, {1,0,2}, {1,2,0}, {2,1,0}.
 */

/* History
 * 08/20/03 (DAK) new
 */
 
public class PermutationIterator {

	private final PermutationIterator subPermutation;
	private boolean hasNext = false;
	private int[] subIndex;
	private int iLast;//place to position largest element
	private int n;
	
	//constructor for use by NullSinglet inner class
	private PermutationIterator() {
		subPermutation = null;
	}
	public PermutationIterator(int n) {
		super();
		if(n < 0) throw new IllegalArgumentException("PermutationIterator constructor requires positive argument");
		this.n = n;
		subPermutation = (n==1) ? new NullSinglet() : new PermutationIterator(n-1);
		reset();
	}

	public void reset() {
		subPermutation.reset();
		subIndex = subPermutation.next();
		iLast = n-1;
		hasNext = true;
	}
	
	public boolean hasNext() {
		return hasNext;
	}
	
	public int[] next() {
		int[] nextIndex = new int[n];
		int j = 0;
		for(int i=0; i<n; i++) {
			if(i==iLast) continue;
			else nextIndex[i] = subIndex[j++];
		}
		nextIndex[iLast--] = n-1;
		if(iLast < 0) {
			if(subPermutation.hasNext()) {
				subIndex = subPermutation.next();
				iLast = n-1;
			} else {
				hasNext = false;
			}
		}
		return nextIndex;
	}
	
	/**
	 * Returns a new matrix formed by permuting the rows and columns of the
	 * given matrix according to the permuation list p.  The list p[i] is such
	 * that it gives the (row/column) index in matrix a from which the element
	 * of the new matrix is taken.  So b(p[i],p[j]) will be equal to a(i,j). The
	 * original matrix a is unchanged.
	 * @param p permutation list
	 * @param a original matrix
	 * @return int[][] permuted matrix
	 */
	public static int[][] matrixPermutation(int[] p, int[][] a) {
		int n = p.length;
		if(a.length != n) throw new IllegalArgumentException("Error: permutation sequence and matrix for permutation are not of compatible size");
		int[][] b = new int[n][n];
		for(int i=0; i<n; i++) {
			int ip = p[i];
			for(int j=0; j<n; j++) {
				b[ip][p[j]] = a[i][j];
			}
		}
		return b;
	}
		
	// Iterator that returns one null value (which isn't used) then expires.  Used to close recursive
	// implemention of PermutationIterator.
	private static class NullSinglet extends PermutationIterator {
		private boolean hasNext = false;
		public boolean hasNext() {return hasNext;}
		public void reset() {hasNext = true;}
		public int[] next() {
			hasNext = false;
			return null;
		}
	}
	
	public static void main(String[] args) {
		int n = Integer.parseInt(args[0]);
		PermutationIterator p = new PermutationIterator(n);
		while(p.hasNext) {
			int[] a = p.next();
			System.out.print("{");
			for(int i=0; i<a.length-1; i++) System.out.print(a[i]+", ");
			System.out.println(a[a.length-1]+"}");
		}
	}
}
