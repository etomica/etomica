package etomica.virial.cluster;
import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

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
	
	/** Prepares a set of bond pairs for a clusters
	 * formed from a disconnected set of fully connected subclusters.
	 * @param iSet element [k] describes the number of subclusters having k
	 * points
	 * @return int[][] the bond array built to the given specification
	 */  
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
	
	public static Cluster[] B6Clusters(MayerFunction f) {

		int[][] FRH1 = new int[][] {{0,1},{0,2},{0,4},{0,5},{1,2},{1,3},{1,5}
										,{2,3},{2,4},{2,5},{3,4},{3,5},{4,5}};
		int[][] FRH2 = new int[][] {{0,2},{0,3},{0,4},{0,5},{1,2},{1,3},{1,4},{1,5}
										,{2,4},{2,5},{3,4},{3,5}};
		int[][] FRH3 = new int[][] {{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{1,5}
										,{2,3},{2,5},{3,4},{3,5},{4,5}};
		int[][] FRH4 = new int[][] {{0,1},{0,2},{0,3},{0,4},{1,3},{1,4},{1,5}
										,{2,4},{2,5},{3,4},{3,5}};
		int[][] FRH5 = new int[][] {{0,2},{0,3},{0,4},{0,5},{1,2},{1,3},{1,4},{1,5}
										,{2,4},{2,5},{3,5}};										
		int[][] FRH6 = new int[][] {{0,1},{0,2},{0,3},{1,2},{1,4},{1,5}
										,{2,3},{2,4},{2,5},{3,4},{3,5}};
										
		int[][] FRH7 = new int[][] {{0,1},{0,3},{0,4},{0,5},{1,2},{1,3},{1,5}
										,{2,3},{2,4},{2,5},{3,5}};
		int[][] FRH8 = new int[][] {{0,2},{0,3},{0,4},{1,3},{1,4},{1,5}
										,{2,4},{2,5},{3,4},{3,5}};
		int[][] FRH9 = new int[][] {{0,2},{0,3},{0,4},{0,5},{1,2},{1,3},{1,4},{1,5}
										,{2,4},{3,5}};		
		int[][] FRH10 = new int[][] {{0,2},{0,3},{0,4},{1,3},{1,4},{1,5}
										,{2,3},{2,5},{3,4},{3,5}};	
		int[][] FRH11 = new int[][] {{0,1},{0,3},{0,4},{0,5},{1,2},{1,3},{1,5}
										,{2,3},{2,4},{2,5}};	
		int[][] FRH12 = new int[][] {{0,2},{0,3},{0,4},{1,3},{1,4},{1,5}
										,{2,4},{2,5},{3,5}};		
		int[][] FRH13 = new int[][] {{0,2},{0,3},{1,3},{1,4},{1,5}
										,{2,4},{2,5},{3,4},{3,5}};		
		int[][] FRH14 = new int[][] {{0,2},{0,3},{0,4},{1,3},{1,4},{1,5}
										,{2,5},{3,4},{3,5}};	
		int[][] FRH15 = new int[][] {{0,2},{0,3},{0,4},{1,4},{1,5}
										,{2,4},{2,5},{3,5}};	
		int[][] FRH16 = new int[][] {{0,2},{0,3},{1,4},{1,5}
										,{2,4},{2,5},{3,4},{3,5}};	
//error		int[][] FRH17 = new int[][] {{0,2},{0,3},{0,4},{1,3},{1,4},{1,5}
//error										,{2,5},{3,5}};		
		int[][] FRH17 = new int[][] {{0,2},{0,3},{0,4},{1,3},{1,4},{1,5}
										,{2,5},{3,4}};		
		int[][] FRH18 = new int[][] {{0,2},{0,3},{0,4},{1,3},{1,4},{1,5}
										,{2,5}};
		int[][] FRH19 = new int[][] {{0,2},{0,3},{1,4},{1,5}
										,{2,4},{3,5}};		
										
		// rest which are negligible in hard spheres
		int[][] FRH20 = new  int[][]{{0,1},{0,2},{0,3},{1,4},{1,5},{2,4},{2,5}
										,{3,4},{3,5}};
		
		int[][] FRH21  = new int[][]{{0,2},{0,3},{0,4},{0,5},{1,2},{1,3},{1,4},{1,5}
										,{2,4}};
										
//error		int[][] FRH22  = new int[][]{{0,2},{0,3},{0,4},{0,5},{1,3},{1,4},{1,5}};						
		int[][] FRH22  = new int[][]{{0,2},{0,3},{0,4},{0,5},{1,2},{1,3},{1,4},{1,5}};						
		
										
																																																														
		Cluster f0 = new Cluster(6, -24.0/144.0, new Cluster.BondGroup(f, Standard.full(6)));
		
		Cluster f1 = new ReeHoover(6, 540.0/144.0, new Cluster.BondGroup(f, FRH1));
		Cluster f2 = new ReeHoover(6, -240/144.0, new Cluster.BondGroup(f, FRH2));
		Cluster f3 = new ReeHoover(6, -1440.0/144.0, new Cluster.BondGroup(f, FRH3));
		Cluster f4 = new ReeHoover(6, 360.0/144.0, new Cluster.BondGroup(f, FRH4));
		Cluster f5 = new ReeHoover(6, 900.0/144.0, new Cluster.BondGroup(f, FRH5));
		Cluster f6 = new ReeHoover(6, 240.0/144.0, new Cluster.BondGroup(f, FRH6));
		Cluster f7 = new ReeHoover(6, 360.0/144.0, new Cluster.BondGroup(f, FRH7));
		Cluster f8 = new ReeHoover(6, 360.0/144.0, new Cluster.BondGroup(f, FRH8));
		Cluster f9 = new ReeHoover(6, -180.0/144.0, new Cluster.BondGroup(f, FRH9));	
		Cluster f10 = new ReeHoover(6, 288.0/144.0, new Cluster.BondGroup(f, FRH10));	
		Cluster f11 = new ReeHoover(6, -540.0/144.0, new Cluster.BondGroup(f, FRH11));		
		Cluster f12 = new ReeHoover(6, -240.0/144.0, new Cluster.BondGroup(f, FRH12));	
		Cluster f13 = new ReeHoover(6, -360.0/144.0, new Cluster.BondGroup(f, FRH13));	
		Cluster f14 = new ReeHoover(6, -1080.0/144.0, new Cluster.BondGroup(f, FRH14));
		Cluster f15 = new ReeHoover(6, 720.0/144.0, new Cluster.BondGroup(f, FRH15));	
		Cluster f16 = new ReeHoover(6, 90.0/144.0, new Cluster.BondGroup(f, FRH16));		
		Cluster f17 = new ReeHoover(6, 360.0/144.0, new Cluster.BondGroup(f, FRH17));	
		Cluster f18 = new ReeHoover(6, -180.0/144.0, new Cluster.BondGroup(f, FRH18));		
		Cluster f19 = new ReeHoover(6, -60.0/144.0, new Cluster.BondGroup(f, FRH19));	
		Cluster f20 = new ReeHoover(6, -40.0/144.0, new Cluster.BondGroup(f, FRH20));	
		Cluster f21 = new ReeHoover(6, 180.0/144.0, new Cluster.BondGroup(f, FRH21));	
		Cluster f22 = new ReeHoover(6, -15.0/144.0, new Cluster.BondGroup(f, FRH22));	
				
		Cluster[] clusters = new Cluster[] {f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22};
		for(int i=0; i<clusters.length; i++) clusters[i].setUsePermutations(true);
		return clusters;
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
