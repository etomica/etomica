package etomica.virial.cluster;

/**
 * @author kofke
 *
 * Class that provides some standard pair sets (via static methods or fields of
 * integer arrays) used in specification of clusters.
 */
public final class Standard {

	/**
	 * Private constructor to prevent instantiation.
	 */
	private Standard() {
		super();
	}


	/**
	 * Returns a chain of bonds, {{0,1},{1,2},...{n-2,n-1}}
	 * @param n number of points in chain
	 * @return int[][] array describing chain of bonds
	 */
	public static int[][] chain(int n) {
		int[][] array = new int[n-1][];
		for(int i=0; i<n-1; i++) {
			array[i] = new int[] {i,i+1};
		}
		return array;
	}
	
	
	/**
	 * Returns a ring of bonds, {{0,1},{1,2},...{n-2,n-1},{n-1,0}}
	 * @param n number of points in ring
	 * @return int[][] array describing ring of bonds
	 */
	public static int[][] ring(int n) {
		int[][] array = new int[n][];
		for(int i=0; i<n-1; i++) {
			array[i] = new int[] {i,i+1};
		}
		array[n-1] = new int[] {0,n-1};
		return array;
	}
	
	/**
	 * Returns a full set of bonds, such that each of the n points is joined to
	 * each of the others.  Starts labeling points at 0.
	 */
	public static int[][] full(int n) {
		return full(n, 0);
	}
	
	/**
	 * Returns a full set of bonds, such that each of the n points is joined to
	 * each of the others; starts labeling of point using the given index
	 * <first>.  For example, full(3,2) returns {{2,3},{2,4},{3,4}}, which can
	 * be compared to full(3), which returns {{0,1},{0,2},{1,2}}
	 */
	public static int[][] full(int n, int first) {
		int[][] array = new int[n*(n-1)/2][];
		int k = 0;
		for(int i=0; i<n-1; i++) {
			for(int j=i+1; j<n; j++) {
				array[k++] = new int[] {i+first,j+first};
			}
		}
		return array;
	}
	
	public static int[][] product(int[] iSet) {
		int n = 0;
		for(int i=1; i<=iSet.length; i++) n += iSet[i-1]*i*(i-1)/2;
		int[][] array = new int[n][];
		int j = 0;
		int first = iSet[0];//only single points (Q1) (no integer pairs) for number of points equal to iSet[0]
		for(int i=2; i<=iSet.length; i++) { //{1, 2, 0, 0, 0} = {Q1, 2Q2, etc}  start with Q2
			for(int k=0; k<iSet[i-1]; k++) {//e.g. make 2 Q2 sets {0,1}, {2,3}
				int[][] a = full(i, first);
				for(int m=0; m<a.length; m++) array[j++] = a[m];//add integer pairs for this group to array of all pairs
				first += i;
			}
		}
		return array;
	}

	public static final int[][] B2 = new int[][] {{0,1}};
	public static final int[][] C3 = ring(3);
	
	public static final int[][] D4 = ring(4);
	public static final int[][] D5 = new int[][] {{0,1},{0,2},{0,3},{1,2},{2,3}};
	public static final int[][] D6 = full(4);
	
	public static double B2HS(double sigma) {
		return 2.0*Math.PI/3.0 * sigma*sigma*sigma;
	}
	public static double C3HS(double sigma) {
		double b0 = B2HS(sigma);
		return 5./8. * b0 * b0;
	}
	
	public static void main(String[] args) {
		int[] iSet = new int[] {5,0,0,0,0};
		int[][] array = product(iSet);		
		int n = array.length;
			String string = "(";
			for(int i=0; i<n-1; i++) string += "{"+array[i][0]+","+array[i][1]+"}, ";
			if(n>0) string += "{"+array[n-1][0]+","+array[n-1][1]+"}";
			string += ")";
		System.out.println(string);
	}
	
//	24 (5, 0, 0, 0, 0)
//	-60 (3, 1, 0, 0, 0)
//	30 (1, 2, 0, 0, 0)
//	20 (2, 0, 1, 0, 0)
//	-10 (0, 1, 1, 0, 0)
//	-5 (1, 0, 0, 1, 0)
//	1 (0, 0, 0, 0, 1)
}
