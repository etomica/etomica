package etomica.virial.cluster;
import etomica.math.SpecialFunctions;
import etomica.util.Arrays;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumEF;
import etomica.virial.ClusterSumStickyEF;
import etomica.virial.ClusterTree;
import etomica.virial.ClusterTreeEF;
import etomica.virial.MayerFunction;
import etomica.virial.MayerHardSphere;

/**
 * @author kofke
 *
 * Class that provides some standard pair sets (via static methods or fields of
 * integer arrays) used in specification of clusters.
 */
public final class Standard implements java.io.Serializable {

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

    public static ClusterAbstract virialCluster(int nBody, MayerFunction f) {
        return virialCluster(nBody,f,true,new FTilde(f),true, true);
    }
    public static ClusterAbstract virialCluster(int nBody, MayerFunction f, 
            boolean usePermutations, MayerFunction e, boolean useTree, boolean uniqueOnly) {
        uniqueOnly = uniqueOnly && nBody > 3 && !useTree;
        if (nBody < 4) {
            e = null;
        }
        if (useTree && e==null) {
            useTree=false;
        }
        if (useTree) {
            System.out.println("using cluster Tree");
            //permutations don't make sense with a tree
            usePermutations = false;
        }
        ClusterDiagramTree tree = null;
        if (useTree) {
            tree = new ClusterDiagramTree(nBody,2);
        }
        int nBondTypes = (e == null) ? 1 : 2;
        ClusterDiagram clusterD = new ClusterDiagram(nBody,0);
        ClusterGenerator generator = new ClusterGenerator(clusterD);
        generator.setAllPermutations(!usePermutations);
        generator.setOnlyDoublyConnected(true);
        generator.setExcludeArticulationPoint(false);
        generator.setExcludeArticulationPair(false);
        generator.setExcludeNodalPoint(false);
        generator.setMakeReeHover(e != null);
        clusterD.reset();
        generator.reset();
        if (e != null) {
            generator.calcReeHoover();
        }
        ClusterBonds[] clusters = new ClusterBonds[0];
        double[] weights = new double[0];
        int fullSymmetry = usePermutations ? 1 : SpecialFunctions.factorial(nBody);
        double weightPrefactor = (1-nBody)/(double)fullSymmetry;
        if (useTree) {
            tree.setPrefactor(weightPrefactor);
        }
        do {
            if (useTree) {
                tree.addClusterDiagram(clusterD);
                continue;
            }
            int iBond = 0, iEBond = 0;
            int numBonds = clusterD.getNumConnections();
            int[][][] bondList = new int[nBondTypes][][];
            bondList[0] = new int[numBonds][2];
            if (nBondTypes == 2) {
                int totalBonds = nBody*(nBody-1)/2;
                bondList[1] = new int[totalBonds-numBonds][2];
            }
            for (int i = 0; i < nBody; i++) {
                int lastBond = i;
                int[] iConnections = clusterD.mConnections[i];
                for (int j=0; j<nBody-1; j++) {
                    if (iConnections[j] > i) {
                        if (e != null) {
                            
                            for (int k=lastBond+1; k<iConnections[j]; k++) {
                                bondList[1][iEBond][0] = i;
                                bondList[1][iEBond++][1] = k;
                            }
                        }
                        bondList[0][iBond][0] = i;
                        bondList[0][iBond++][1] = iConnections[j];
                        lastBond = iConnections[j];
                    }
                    else if ((lastBond>i || iConnections[j] == -1) && e != null) {
                        for (int k=lastBond+1; k<nBody; k++) {
                            bondList[1][iEBond][0] = i;
                            bondList[1][iEBond++][1] = k;
                        }
                    }
                    if (iConnections[j] == -1) break;
                }
            }
            // only use permutations if the diagram has permutations
            boolean thisUsePermutations = !uniqueOnly && usePermutations && 
                                          clusterD.mNumIdenticalPermutations < SpecialFunctions.factorial(nBody);
            // only use e-bonds if one of the diagrms has some
            clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, thisUsePermutations));
            double [] newWeights = new double[weights.length+1];
            System.arraycopy(weights,0,newWeights,0,weights.length);
            newWeights[weights.length] = clusterD.mReeHooverFactor*weightPrefactor/clusterD.mNumIdenticalPermutations;
            weights = newWeights;
        } while (generator.advance());
        if (true) {
            if (useTree) {
                tree.collapse();
                return new ClusterTreeEF(tree,new MayerFunction[]{e});
            }
            if (e != null) {
                return new ClusterSumEF(clusters,weights,new MayerFunction[]{e});
            }
            return new ClusterSum(clusters,weights,new MayerFunction[]{f});
        }
        MayerFunction[] allF = new MayerFunction[nBondTypes];
        allF[0] = f;
        if (e != null) {
            allF[1] = e;
        }
        if (useTree) {
            tree.collapse();
            return new ClusterTree(tree,allF);
        }
        return new ClusterSum(clusters,weights,allF);
    }
    
	public static final int[][] B2 = new int[][] {{0,1}};
	public static final int[][] C3 = ring(3);
	
	public static final int[][] D4 = ring(4);
	public static final int[][] D5 = new int[][] {{0,1},{0,2},{0,3},{1,2},{2,3}};
	public static final int[][] D6 = full(4);
	
	public static double B2HS(double sigma) {
		return 2.0*Math.PI/3.0 * sigma*sigma*sigma;
	}
	public static double B3HS(double sigma) {
		double b0 = B2HS(sigma);
		return 5./8. * b0 * b0;
	}
    
    public static double B4HS(double sigma) {
        double b0 = B2HS(sigma);
        return (219.0*Math.sqrt(2.0)/2240.0/Math.PI-89.0/280.0+4131.0/2240.0/Math.PI*Math.atan(Math.sqrt(2.0)))*b0*b0*b0;
    }
	
    public static double B5HS(double sigma) {
        double b0 = B2HS(sigma);
        return 0.110252*b0*b0*b0*b0;
    }
    
    public static double B6HS(double sigma) {
        double b0 = B2HS(sigma);
        return 0.038808*b0*b0*b0*b0*b0;
    }

    public static double B7HS(double sigma) {
        double b0 = B2HS(sigma);
        return 0.013046*b0*b0*b0*b0*b0*b0;
    }
    
    public static double B2SW(double sigma, double lambda, double ekT) {
        if (lambda < 1) {
            return B2HS(sigma);
        }
        double d = Math.exp(ekT)-1;
        return B2HS(sigma)*(1-(lambda*lambda*lambda-1)*d);
    }

    public static double B3SW(double sigma, double lambda, double ekT) {
        if (lambda < 1) {
            return B3HS(sigma);
        }
        double d = Math.exp(ekT)-1;
        double f1, f2, f3;
        double lambda2 = lambda*lambda;
        double lambda3 = lambda*lambda*lambda;
        double lambda4 = lambda2*lambda2;
        double lambda6 = lambda3*lambda3;
        if (lambda < 2) {
            f1 = 0.2*(lambda6 - 18*lambda4 + 32*lambda3 - 15);
            f2 = 0.4*(lambda6 - 18*lambda4 + 16*lambda3 + 9*lambda2 - 8);
            f3 = 1.2*Math.pow(lambda2-1,3);
        }
        else {
            f1 = 3.4;
            f2 = 0.2*(          - 32*lambda3 + 18*lambda2 + 48);
            f3 = 0.2*(5*lambda6 - 32*lambda3 + 18*lambda2 + 26);
        }
        return B3HS(sigma)*(1-f1*d-f2*d*d-f3*d*d*d);
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
}
