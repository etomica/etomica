/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.math.Rational;
import etomica.math.SpecialFunctions;
import etomica.util.Arrays;
import etomica.virial.MayerFunction;

import java.util.ArrayList;

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

    public static ClusterAbstract virialCluster(int nBody, MayerFunction f, 
            boolean usePermutations, MayerFunction e, boolean uniqueOnly) {
        uniqueOnly = uniqueOnly && nBody > 3;
        if (nBody < 4) {
            e = null;
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
        long fullSymmetry = usePermutations ? 1 : SpecialFunctions.factorial(nBody);
        double weightPrefactor = (1-nBody)/(double)fullSymmetry;
        do {
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
            // only use e-bonds if one of the diagram has some
            clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, thisUsePermutations));
            double [] newWeights = new double[weights.length+1];
            System.arraycopy(weights,0,newWeights,0,weights.length);
            newWeights[weights.length] = clusterD.mReeHooverFactor*weightPrefactor/clusterD.mNumIdenticalPermutations;
            weights = newWeights;
        } while (generator.advance());
        return new ClusterSum(clusters,weights,new MayerFunction[]{f});
    }
    
    public static ClusterAbstract virialSeriesCluster(int nBody, MayerFunction[] f, int[] bondTypes) {
        ClusterDiagram clusterD = new ClusterDiagram(nBody,0);
        ClusterGenerator generator = new ClusterGenerator(clusterD);
        generator.setAllPermutations(false);
        generator.setOnlyDoublyConnected(true);
        generator.setExcludeArticulationPoint(false);
        generator.setExcludeArticulationPair(false);
        generator.setExcludeNodalPoint(false);
        generator.setMakeReeHover(false);
        clusterD.reset();
        generator.reset();
        int nBondTypes = f.length;
        PermutationIterator iter = new PermutationIterator(bondTypes);
        int numBonds = bondTypes.length;
        ClusterBonds[] clusters = new ClusterBonds[0];
        double[] weights = new double[0];
        while (true) {
            while (clusterD.getNumConnections() != numBonds) {
                if (!generator.advance()) {
                    if (clusters.length > 0) {
                        return new ClusterSum(clusters, weights, f);
                    }
                    throw new RuntimeException("couldn't find number of bonds "+bondTypes.length);
                }
            }
            iter.reset();
            int nPermutations = 0;
            for (int[] iBondTypes = iter.next(); iBondTypes != null; iBondTypes = iter.next()) {
                nPermutations++;
            }
            iter.reset();
            for (int[] iBondTypes = iter.next(); iBondTypes != null; iBondTypes = iter.next()) {
                int[][][] bondList = new int[nBondTypes][][];
                for (int i=0; i<nBondTypes; i++) {
                    int n = 0;
                    for (int j=0; j<numBonds; j++) {
                        if (iBondTypes[j] == i) {
                            n++;
                        }
                    }
                    bondList[i] = new int[n][2];
                }
                int bondCount = 0;
                int[] iBondCount = new int[nBondTypes];
                for (int i = 0; i < nBody; i++) {
                    int[] iConnections = clusterD.mConnections[i];
                    for (int j=0; j<nBody-1 && iConnections[j] != -1; j++) {
                        if (iConnections[j] > i) {
                            int bondType = iBondTypes[bondCount];
                            bondList[bondType][iBondCount[bondType]][0] = i;
                            bondList[bondType][iBondCount[bondType]][1] = iConnections[j];
                            iBondCount[bondType]++;
                            bondCount++;
                        }
                    }
                }
                clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
                double [] newWeights = new double[weights.length+1];
                System.arraycopy(weights,0,newWeights,0,weights.length);
                newWeights[weights.length] = (1.0-nBody)/clusterD.mNumIdenticalPermutations/nPermutations;
                System.out.println((1.0-nBody)+" "+clusterD.mNumIdenticalPermutations);
                weights = newWeights;
            }
            if (!generator.advance()) {
                if (clusters.length > 0) {
                    return new ClusterSum(clusters, weights, f);
                }
                throw new RuntimeException("couldn't find number of bonds "+bondTypes.length);
            }
        }
    }
    
    public static ClusterAbstract virialClusterXS(int nBody, MayerFunction f, 
            boolean usePermutations, MayerFunction e, boolean uniqueOnly, int approx) {
        uniqueOnly = uniqueOnly && nBody > 3;
        if (nBody < 4) {
            e = null;
        }
        int nBondTypes = (e == null) ? 1 : 2;
        ClusterDiagram clusterD = new ClusterDiagram(nBody,2);
        ClusterGenerator generator = new ClusterGenerator(clusterD);
        generator.setAllPermutations(false);
        generator.setOnlyConnected(false);
        generator.setOnlyDoublyConnected(true);
        generator.setExcludeArticulationPoint(true);
        generator.setExcludeArticulationPair(false);
        generator.setExcludeNodalPoint(true);
        generator.setMakeReeHover(false);
//        clusterD.reset();
        generator.reset();


        clusterD.setWeight(new Rational(1, clusterD.mNumIdenticalPermutations));
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        list.add(new ClusterDiagram(clusterD));
        while(generator.advance()) {
            clusterD.setWeight(new Rational(1, clusterD.mNumIdenticalPermutations));
            list.add(new ClusterDiagram(clusterD));
        }
        ClusterDiagram[] trueClusters = list.toArray(new ClusterDiagram[list.size()]);
        System.out.println("true clusters");
        for (int i=0; i<trueClusters.length; i++) {
            System.out.println(trueClusters[i]);
        }
        ClusterOperations.addEquivalents(list);
        trueClusters = list.toArray(new ClusterDiagram[list.size()]);
        System.out.println("true clusters (equiv)");
        for (int i=0; i<trueClusters.length; i++) {
            System.out.println(trueClusters[i]);
        }
        ClusterDiagram[] out = trueClusters;
        if (approx != ClusterOperations.NONE) {
            ClusterOperations ops = new ClusterOperations();
            ops.setApproximation(approx);
            ClusterDiagram[] approxClusters = ops.getC(nBody-2);
            System.out.println("approx clusters");
            for (int i=0; i<approxClusters.length; i++) {
                System.out.println(approxClusters[i]);
            }
            out = ClusterOperations.difference(trueClusters, approxClusters);
            System.out.println("difference clusters");
            for (int i=0; i<out.length; i++) {
                System.out.println(out[i]);
            }
        }
        out = ClusterOperations.integrate(out);
        System.out.println("integrate1 clusters");
        for (int i=0; i<out.length; i++) {
            System.out.println(out[i]);
        }
        out = ClusterOperations.integrate(out);
        System.out.println("integrate2 clusters");
        for (int i=0; i<out.length; i++) {
            System.out.println(out[i]);
        }
        if (e != null) {
            out = ClusterOperations.makeReeHoover(out);
            System.out.println("ree-hoovered clusters");
            for (int i=0; i<out.length; i++) {
                System.out.println(out[i]);
            }
        }
//        for (int i=0; i<xs.length; i++) {
//            System.out.println(xs[i].toString());
//        }
//        System.exit(1);

        ClusterBonds[] clusters = new ClusterBonds[0];
        double[] weights = new double[0];
        long fullSymmetry = usePermutations ? 1 : SpecialFunctions.factorial(nBody);
        double weightPrefactor = (1-nBody)/(double)fullSymmetry;

        for (int m=0; m<out.length; m++) {
            clusterD = out[m];
            ClusterOperations.sortConnections(clusterD);
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
                        // we're done with f-bonds, fill in the rest with e-bonds
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
            // only use e-bonds if one of the diagrams has some
            clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, thisUsePermutations));
            double [] newWeights = new double[weights.length+1];
            System.arraycopy(weights,0,newWeights,0,weights.length);
            newWeights[weights.length] = clusterD.mReeHooverFactor*weightPrefactor/clusterD.mNumIdenticalPermutations;
            System.out.println("hi "+clusterD.mReeHooverFactor+" "+weightPrefactor+" "+clusterD.mNumIdenticalPermutations+" "+newWeights[weights.length]);
            weights = newWeights;
        } while (generator.advance());
        System.out.println("XS weights: "+java.util.Arrays.toString(weights));
        return new ClusterSum(clusters,weights,new MayerFunction[]{f});
    }
    
    public static ClusterAbstract virialClusterMixture(int nBody, MayerFunction[][] f, MayerFunction[][] e, int[] nTypes) {
        if (nBody < 4) {
            e = null;
        }
        int[] pointType = new int[nBody];
        int l = 0;
        // label the first points to be type 1, the next points to be type 2, etc
        for (int i=0; i<nTypes.length; i++) {
            for (int j=0; j<nTypes[i]; j++) {
                pointType[l] = i;
                l++;
            }
        }
        // nTypes.length is the number of components
        // we need one bond type for each pair of components
        int nBondTypes = nTypes.length*(nTypes.length+1)/2;
        // bondType is bond index for the type of bond between points of type i and j
        int[][] bondType = new int[nTypes.length][nTypes.length];
        // a linear list of the f and e functions (which come in as 2D, indexed by i and j)
        MayerFunction[] linearF = new MayerFunction[nBondTypes];
        MayerFunction[] linearE = null;
        if (e != null) {
            linearE = new MayerFunction[nBondTypes];
        }
        l=0;
        for (int i=0; i<nTypes.length; i++) {
            for (int j=0; j<i+1; j++) {
                bondType[i][j] = l;
                // we're symmetric
                bondType[j][i] = l;
                linearF[l] = f[i][j];
                if (e != null) {
                    linearE[l] = e[i][j];
                }
                // we ignore f[j][i] and e[j][i] since e and f should be symmetric.
                l++;
            }
        }
        int allNumBondTypes = nBondTypes;
        if (e != null) {
            // with e-bonds, we need an extra bond type for each pair of components
            allNumBondTypes *= 2;
        }
        ClusterDiagram clusterD = new ClusterDiagram(nBody,0);
        ClusterGenerator generator = new ClusterGenerator(clusterD);
        generator.setAllPermutations(true);
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
        long fullSymmetry = SpecialFunctions.factorial(nBody);
        double weightPrefactor = (1-nBody)/(double)fullSymmetry;
        int[] iBond = new int[allNumBondTypes];
        do {
            for (int i=0; i<allNumBondTypes; i++) {
                iBond[i] = 0;
            }
            int numBonds = clusterD.getNumConnections();
            // bondList[i][j][0] is the first point for the jth bond of type i
            // bondList[i][j][1] is the second point for the jth bond of type i
            int[][][] bondList = new int[allNumBondTypes][][];
            for (l=0; l<allNumBondTypes; l++) {
                bondList[l] = new int[numBonds][2];
            }
            if (e != null) {
                int totalBonds = nBody*(nBody-1)/2;
                for (l=0; l<nBondTypes; l++) {
                    bondList[nBondTypes+l] = new int[totalBonds-numBonds][2];
                }
            }
            for (int i = 0; i < nBody; i++) {
                int lastBond = i;
                int[] iConnections = clusterD.mConnections[i];
                for (int j=0; j<nBody-1; j++) {
                    if (iConnections[j] > i) {
                        if (e != null) {
                            
                            for (int k=lastBond+1; k<iConnections[j]; k++) {
                                // thisBondType is the bond type connecting point i and k
                                int thisBondType = nBondTypes+bondType[pointType[i]][pointType[k]];
                                bondList[thisBondType][iBond[thisBondType]][0] = i;
                                bondList[thisBondType][iBond[thisBondType]++][1] = k;
                            }
                        }
                        
                        // thisBondType is the bond type connecting point i and iConnections[j]
                        int thisBondType = bondType[pointType[i]][pointType[iConnections[j]]];
                        bondList[thisBondType][iBond[thisBondType]][0] = i;
                        bondList[thisBondType][iBond[thisBondType]++][1] = iConnections[j];
                        lastBond = iConnections[j];
                    }
                    else if ((lastBond>i || iConnections[j] == -1) && e != null) {
                        for (int k=lastBond+1; k<nBody; k++) {
                            int thisBondType = nBondTypes+bondType[pointType[i]][pointType[k]];
                            bondList[thisBondType][iBond[thisBondType]][0] = i;
                            bondList[thisBondType][iBond[thisBondType]++][1] = k;
                        }
                    }
                    if (iConnections[j] == -1) break;
                }
            }
            // we oversized bondList because we didn't know how much we'd need.
            // now resize the bondList for each bondType
            for (int i=0; i<bondList.length; i++) {
                if (iBond[i] == 0) {
                    bondList[i] = new int[0][0];
                    continue;
                }
                int[][] newBondList = new int[iBond[i]][2];
                System.arraycopy(bondList[i], 0, newBondList, 0, iBond[i]);
                bondList[i] = newBondList;
            }
            clusters = (ClusterBonds[])Arrays.addObject(clusters,new ClusterBonds(nBody, bondList, false));
            double [] newWeights = new double[weights.length+1];
            System.arraycopy(weights,0,newWeights,0,weights.length);
            newWeights[weights.length] = clusterD.mReeHooverFactor*weightPrefactor/clusterD.mNumIdenticalPermutations;
            weights = newWeights;
        } while (generator.advance());
        return new ClusterSum(clusters,weights,linearF);
    }

	public static final int[][] B2 = new int[][] {{0,1}};
	public static final int[][] C3 = ring(3);
	
	public static final int[][] D4 = ring(4);
	public static final int[][] D5 = new int[][] {{0,1},{0,2},{0,3},{1,2},{2,3}};
	public static final int[][] D6 = full(4);
	
	/**
	 * Returns the nth virial coefficient for hard spheres of diameter sigma.
	 */
	public static double BHS(int n, double sigma) {
	    switch (n) {
	        case 2:
	            return B2HS(sigma);
	        case 3:
	            return B3HS(sigma);
            case 4:
                return B4HS(sigma);
            case 5:
                return B5HS(sigma);
            case 6:
                return B6HS(sigma);
            case 7:
                return B7HS(sigma);
            case 8:
                return B8HS(sigma);
            case 9:
                return B9HS(sigma);
            case 10:
                return B10HS(sigma);
            case 11:
                return B11HS(sigma);
            case 12:
                return B12HS(sigma);
            default:
                throw new RuntimeException("HS unknown at order "+n);
	    }
	}
	
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
        // Labik, Kolafa and Malijevsky, Phys. Rev. E (71), 2005, 021105
        // 28.22445(10)
        return 28.22445*Math.pow(b0/4, 4);
    }
    
    public static double B6HS(double sigma) {
        double b0 = B2HS(sigma);
        // Labik, Kolafa and Malijevsky, Phys. Rev. E (71), 2005, 021105
        // 39.81550(36)
        return 39.81550*Math.pow(b0/4, 5);
    }

    public static double B7HS(double sigma) {
        double b0 = B2HS(sigma);
        // Labik, Kolafa and Malijevsky, Phys. Rev. E (71), 2005, 021105
        // 53.3413(16)
        return 53.3413*Math.pow(b0/4, 6);
    }
    
    public static double B8HS(double sigma) {
        double b0 = B2HS(sigma);
        // Labik, Kolafa and Malijevsky, Phys. Rev. E (71), 2005, 021105
        // 68.540(10) / 4^7
        return 68.540*Math.pow(b0/4, 7);
    }
    
    public static double B9HS(double sigma) {
        double b0 = B2HS(sigma);
        // Clisby and McCoy, J. Statistical Physics (122) 2006, 15
        // 0.0013094(13)
        return 0.0013094*Math.pow(b0,8);
    }
    
    public static double B10HS(double sigma) {
        double b0 = B2HS(sigma);
        // Clisby and McCoy, J. Statistical Physics (122) 2006, 15
        // 0.0004035(15)
        return 0.0004035*Math.pow(b0,9);
    }
    
    public static double B11HS(double sigma) {
        // Wheatley
        return 0.198*Math.pow(sigma, 30);
    }

    public static double B12HS(double sigma) {
        // Wheatley
        return 0.090*Math.pow(sigma, 33);
    }

    /**
     * Return integral for ring of hard spheres of size nPoints
     */
    public static double ringHS(int nPoints) {
        double[] value = new double[]{0,0,4.0*Math.PI/3.0,
                8.2246703342411321824,    // (5 \[Pi]^2)/6
                23.798821183891502188,    // (2176 \[Pi]^3)/2835
                73.280512726913686846,    // (40949 \[Pi]^4)/54432
                238.01417308225752423,    // (23648512 \[Pi]^5)/30405375
                801.84930772092072806,    // (25040879363 \[Pi]^6)/30023136000
                2776.7662132440473589,    // (35836384927744 \[Pi]^7)/38979295480125
                9823.4092155784055129,    // (6060781370812815989 \[Pi]^8)/5854170457175040000
                35350.421305790123166,    // (430588656377655296 \[Pi]^9)/363092137397364375
                128997.05374244207554,    // (279634343277877995897619523 \[Pi]^10)/203006266387416011520000000
                476213.79444250996002};   // (2181664367287144834983919616 \[Pi]^11)/1347828286825972065254765625
        return value[nPoints];
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
		
		System.out.println(B4HS(1));
		System.out.println(B5HS(1));
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
