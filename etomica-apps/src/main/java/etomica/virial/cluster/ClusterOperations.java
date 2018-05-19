/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedList;

import etomica.math.SpecialFunctions;
import etomica.math.discrete.FixedSumIterator;
import etomica.math.Rational;

public class ClusterOperations {

    public ClusterOperations() {
    }
    
    /**
     * Sorts diagrams in order from most bonds to fewest bonds
     */
    public static void sortDiagrams(ClusterDiagram[] diagrams) {
        java.util.Arrays.sort(diagrams, new Comparator<ClusterDiagram>() {
        
            public int compare(ClusterDiagram o1, ClusterDiagram o2) {
                int n1 = o1.getNumConnections();
                int n2 = o2.getNumConnections();
                return n1 > n2 ? -1 : (n1 < n2 ? 1 : 0);
            }
        });
    }
    
    public static void sortConnections(ClusterDiagram cluster) {
        for (int i=0; i<cluster.mConnections.length; i++) {
            int lastBond = -1;
            for (int j=0; j<cluster.mConnections[i].length; j++) {
                if (cluster.mConnections[i][j] > -1) {
                    lastBond = j;
                }
            }
            java.util.Arrays.sort(cluster.mConnections[i], 0, lastBond+1);
        }
    }
    
    /**
     * Returns a new cluster that is the product of the given clusters.  Input clusters must have the
     * same number of root points and my not both have bonds between the same root points.
     * The weight of the returned cluster is the product of the weights of the given clusters.  The
     * input clusters are not changed.
     */
    public static ClusterDiagram product(ClusterDiagram cluster1, ClusterDiagram cluster2) {
        int numRoot = cluster1.getNumRootPoints();
        if(numRoot != cluster2.getNumRootPoints()) {
            throw new IllegalArgumentException("Product requires diagrams with same number of root points");
        }
        int offset = cluster1.mNumBody - numRoot;
        int numBody = cluster2.mNumBody + offset;
        ClusterDiagram productCluster = new ClusterDiagram(numBody, numRoot, new int[][] {});

        //copy cluster1 bonds to product
        cluster1.copyTo(productCluster);

        //copy cluster2 bonds to product
        for (int i = 0; i < cluster2.mNumBody; i++) {
            int[] connections = cluster2.mConnections[i];
            int node1 = i;
            if(node1 >= numRoot) node1 += offset;
            for(int j = 0; ; j++) {
                int node2 = connections[j];
                if(node2 == -1) {
                    break;
                }
                if(node2 >= numRoot) {
                    node2 += offset;
                }
                productCluster.addConnection(node1, node2);
            }
        }
        productCluster.setWeight(cluster1.getWeight().times(cluster2.getWeight()));
        return productCluster; 
    }
    
    /**
     * Returns the product of all the given clusters.  Length of input array must
     * be at least 1, and no null elements are permitted.  The input array and
     * the clusters it contains are unchanged, and the returned cluster is new.
     */
    public static ClusterDiagram product(ClusterDiagram[] clusters) {
        ClusterDiagram productDiagram = clusters[0];
        for(int i=1; i<clusters.length; i++) {
            productDiagram = product(productDiagram, clusters[i]); 
        }
        return productDiagram;
    }
    
    /**
     * Returns an array of clusters formed from taking the product of each of the diagrams
     * of the first set with each of the diagrams of the second set.  The given arrays and
     * the clusters in them are not changed, and all clusters in the returned array are new.
     */
    public static ClusterDiagram[] product(ClusterDiagram[] set1, ClusterDiagram[] set2) {
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        if(set1.length == 0 || set2.length == 0) {
            return new ClusterDiagram[] {};
        }
        for(int i=0; i<set1.length; i++) {
            for(int j=0; j<set2.length; j++) {
                list.add(product(set1[i],set2[j]));
            }
        }
        addEquivalents(list);
        return list.toArray(new ClusterDiagram[] {});
    }
    
    /**
     * Replaces each cluster in the given array with another having the same
     * topology, but with the highest-index root point converted into a field
     * point. The resulting set is collapsed to add any newly equivalent
     * clusters. Weights of each cluster are otherwise unchanged. The given
     * arrays and the clusters in them are not changed, and all clusters in the
     * returned array are new.
     * 
     * @throws IllegalArgumentException
     *             if any input cluster has no root points
     */
    public static ClusterDiagram[] integrate(ClusterDiagram[] clusters) {
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        for(int i=0; i<clusters.length; i++) {
            if(clusters[i].getNumRootPoints() == 0) {
                throw new IllegalArgumentException("Cannot integrate a cluster having no root points");
            }
            ClusterDiagram newCluster = new ClusterDiagram(clusters[i].mNumBody, 
                    clusters[i].getNumRootPoints()-1);
            clusters[i].copyTo(newCluster);
            list.add(newCluster);
        }
        addEquivalents(list);
        return list.toArray(new ClusterDiagram[] {});
    }
    
    /**
     * Replaces all isomorphic clusters in the given list by a single cluster
     * with weight given by their sum.  The resulting list is also sorted by #
     * of connections.  Both the list and the the weights of the clusters it
     * contains may be altered by this method.
     */
    public static void addEquivalents(ArrayList<ClusterDiagram> list) {
        if (list.size() < 2) {
            return;
        }
        reduce(list);
        ClusterDiagram cluster1 = null;
        for (int i=0; i<list.size(); i++) {
            cluster1 = list.get(i);
            for (int j=i+1; j<list.size(); j++) {
                ClusterDiagram cluster2 = list.get(j);
                if (cluster1.getNumConnections() != cluster2.getNumConnections()) {
                    break;
                }
                if(cluster2.getWeight().numerator() == 0) {
                    continue;
                }
                if(cluster1.isIsomorphOf(cluster2)) {
                    int[] score1 = new int[cluster1.mNumBody/2+1];
                    cluster1.calcScore(score1);
                    if(cluster2.scoreGreaterThan(score1)) {
                        cluster2.setWeight(cluster1.getWeight().plus(cluster2.getWeight()));
                        cluster1.setWeight(Rational.ZERO);
                        break;
                    }
                    cluster1.setWeight(cluster1.getWeight().plus(cluster2.getWeight()));
                    cluster2.setWeight(Rational.ZERO);
                }
            }
        }
        //remove clusters having zero weight
        for (int i=0; i<list.size(); i++) {
            ClusterDiagram cluster = list.get(i);
            if(cluster.getWeight().numerator() == 0) {
                list.remove(i);
                i--;
            }
        }
        list.trimToSize();
    }

    /**
     * Sorts cluster diagrams (but # of connections) and combines any diagrams
     * that are identical, resulting in a diagram with a weight equal to the
     * sum of the original weights.
     */
    public static void reduce(ArrayList<ClusterDiagram> list) {
        if (list.size() < 2) {
            return;
        }
        ClusterDiagram[] array = list.toArray(new ClusterDiagram[list.size()]);
        sortDiagrams(array);

        list.clear();
        int numBody = array[0].mNumBody;
        int[][] scores = new int[array.length][numBody/2+1];
        for (int i=0; i<array.length; i++) {
            ClusterDiagram iCluster = array[i];
            iCluster.calcScore(scores[i]);
            list.add(array[i]);
        }
        for (int i=0; i<list.size()-1; i++) {
            ClusterDiagram iCluster = list.get(i);
            for (int j=i+1; j<list.size(); j++) {
                ClusterDiagram jCluster = list.get(j);
                if (iCluster.getNumConnections() != jCluster.getNumConnections()) {
                    break;
                }
                if (java.util.Arrays.equals(scores[i], scores[j])) {
                    iCluster.setWeight(iCluster.getWeight().plus(jCluster.getWeight()));
                    list.remove(j);
                    for (int k=j; k<scores.length-1; k++) {
                        scores[k] = scores[k+1];
                    }
                    j--;
                    if (iCluster.getWeight().numerator() == 0) {
                        list.remove(i);
                        for (int k=i; k<scores.length-1; k++) {
                            scores[k] = scores[k+1];
                        }
                        i--;
                        break;
                    }                        
                }
            }
        }
    }
    
    /**
     * Returns a new cluster that is a convolution of the given cluster, which is obtained by joining a root point
     * from each diagram and turning it into a field point.  Each given cluster must have exactly 2 root points.
     * The weight of the returned cluster is the product of the weights of the given clusters.  The input clusters
     * are unchanged.
     */
    public static ClusterDiagram convolution(ClusterDiagram cluster1, ClusterDiagram cluster2) {
        int numRoot = cluster1.getNumRootPoints();
        if(numRoot != 2 || numRoot != cluster2.getNumRootPoints()) {
            throw new IllegalArgumentException("Convolution requires diagrams with exactly two root points each");
        }
        int numBody = cluster1.mNumBody + cluster2.mNumBody - 1;
        ClusterDiagram convolutionCluster = new ClusterDiagram(numBody, numRoot, new int[][] {});
        for(int i=0; i<cluster1.mNumBody; i++) {
            int[] connections = cluster1.mConnections[i];
            int node1 = i;
            if(node1 == 0) {
                node1 = cluster1.mNumBody;
            }
            for(int j=0; ; j++) {
                int node2 = connections[j];
                if(node2 == -1) break;
                if(node2 == 0) {
                    node2 = cluster1.mNumBody;
                }
                convolutionCluster.addConnection(node1, node2);
            }
        }
        for(int i=0; i<cluster2.mNumBody; i++) {
            int[] connections = cluster2.mConnections[i];
            int node1 = i;
            if(node1 != 0) {
                node1 += cluster1.mNumBody - 1;
            }
            //if(node1 > 1) node1--;
            for(int j=0; ; j++) {
                int node2 = connections[j];
                if(node2 == -1) break;
                if(node2 != 0) {
                    node2 += cluster1.mNumBody - 1;
                }
                //if(node2 > 1) node2--;
                convolutionCluster.addConnection(node1, node2);
            }
        }
        convolutionCluster.setWeight(cluster1.getWeight().times(cluster2.getWeight()));
        return convolutionCluster;
    }
    
    /**
     * Generates the set of diagrams obtained by subtracting the second set from the first set.
     * The returned array contains all diagrams from the first set, with weights modified by subtracting
     * the weights of matching diagrams from the second set (if present), and all unmatched diagrams from
     * the second set with their weights negated.  The returned array and all clusters in it are new
     * instances.
     */
    public static ClusterDiagram[] difference(ClusterDiagram[] set1, ClusterDiagram[] set2) {
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        for(int i=0; i<set1.length; i++) {
            list.add(new ClusterDiagram(set1[i]));
        }
        for(int i=0; i<set2.length; i++) {
            ClusterDiagram cluster2 = new ClusterDiagram(set2[i]);
            cluster2.setWeight(cluster2.getWeight().times(new Rational(-1,1)));
            list.add(cluster2);
        }
        addEquivalents(list);
        return list.toArray(new ClusterDiagram[] {});
    }

    /**
     * Generates the set of diagrams obtained by adding the second set to the first set.
     * The returned array contains (copies of) all diagrams from the first set, with weights modified by adding
     * the weights of matching diagrams from the second set (if present), and all unmatched diagrams from
     * the second set. The returned array and all clusters in it are new instances.
     */
    public static ClusterDiagram[] sum(ClusterDiagram[] set1, ClusterDiagram[] set2) {
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        for(int i=0; i<set1.length; i++) {
            list.add(new ClusterDiagram(set1[i]));
        }
        for(int i=0; i<set2.length; i++) {
            list.add(new ClusterDiagram(set2[i]));
        }
        addEquivalents(list);
        return list.toArray(new ClusterDiagram[] {});
    }

    /**
     * Returns an array of clusters formed from taking the convolution of each of the diagrams
     * of the first set with each of the diagrams of the second set.  The given arrays and the
     * clusters they change are unchanged, and all clusters in the returned array are new.
     */
    public static ClusterDiagram[] convolution(ClusterDiagram[] set1, ClusterDiagram[] set2) {
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        if(set1.length == 0 || set2.length == 0) {
            return new ClusterDiagram[] {};
        }
        for(int i=0; i<set1.length; i++) {
            for(int j=0; j<set2.length; j++) {
                list.add(convolution(set1[i],set2[j]));
            }
        }
        addEquivalents(list);
        return list.toArray(new ClusterDiagram[] {});
    }
    
    /**
     * Creates a set of Ree-Hoover clusters from the given cluster.  The input cluster is assumed
     * to be formed from only f bonds.  The set of returned clusters are full stars formed from
     * f bonds and (implied) e bonds (e = f + 1), and sum to the original cluster. 
     */
    public static ClusterDiagram[] makeReeHoover(ClusterDiagram cluster) {
        int n = cluster.mNumBody;
        int nUnbonded = n*(n-1)/2 - cluster.getNumConnections();
        ClusterDiagram[] rhClusters = new ClusterDiagram[1<<nUnbonded];
        int[][] pairList = new int[nUnbonded][2];
        int k = -1;
        //make a list of unbonded pairs
        for(int i=0; i<n-1; i++) {
            int[] connections = cluster.mConnections[i];
            boolean[] bonded = new boolean[cluster.mNumBody];
            for(int j=0; connections[j]!=-1; j++) {
                bonded[connections[j]] = true;
            }
            for(int j=i+1; j<n; j++) {
                if(!bonded[j]) {
                    pairList[++k][0] = i;
                    pairList[k][1] = j;
                }
            }
        }
        //bits of i are used to key whether an f (bit = 1) or an e (bit = 0) bond is placed between unbonded pair
        //k is used to index the bonds
        for(int i=0; i<rhClusters.length; i++) {
            rhClusters[i] = new ClusterDiagram(cluster);
            boolean neg = false;
            for(k=0; k<nUnbonded; k++) {
                if((i & (1<<k)) != 0) {//bit for bond k is 1 -- add f-bond
                    neg = !neg;
                    int node1 = pairList[k][0];
                    int node2 = pairList[k][1];
                    rhClusters[i].addConnection(node1, node2);
                    rhClusters[i].addConnection(node2, node1);
                }
            }
            if(neg) {//odd number of f-bonds introduced
                rhClusters[i].setWeight(rhClusters[i].getWeight().times(new Rational(-1,1)));
            }
        }
        return rhClusters;
    }
    
    public static ClusterDiagram[] makeReeHoover(ClusterDiagram[] clusters) {
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        for(int i=0; i<clusters.length; i++) {
            ClusterDiagram[] rhClusters = makeReeHoover(clusters[i]);
            for (int j=0; j<rhClusters.length; j++) {
                if (rhClusters[j].getWeight().numerator() != 0) {
                    list.add(rhClusters[j]);
                }
            }
        }
        addEquivalents(list);
        return list.toArray(new ClusterDiagram[list.size()]);
    }
    
    public void setApproximation(int approximation) {
        approx = approximation;
    }
    
    public int getApproximation() {
        return approx;
    }
    
    /**
     * Returns all diagrams of order n in the expansion for the total correlation function h(r), in the presence
     * of any approximations currently in effect.
     */
    public ClusterDiagram[] getH(int n) {
        if(n < h.length && h[n] != null) {
            return h[n];
        }
        if(n >= h.length) {
            h = Arrays.copyOf(h, n + 1);
        }
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        if(approx == PY) {
            getEta(n);
            for(int i=0; i<eta[n].length; i++) {
                list.add(new ClusterDiagram(eta[n][i]));
                ClusterDiagram eeta = new ClusterDiagram(eta[n][i]);
                eeta.addConnection(0, 1);
                eeta.addConnection(1, 0);
                list.add(eeta);
            }
        } else {
            FixedSumIterator iterator = new FixedSumIterator(n);
            iterator.setSum(n);
            iterator.reset();
            for(int[] arr=iterator.next(); arr!=null; arr=iterator.next()) {
                ClusterDiagram[] clusters = new ClusterDiagram[] {unity};
                int factor = 1;
                for(int i=0; i<n; i++) {
                    factor *= SpecialFunctions.factorial(arr[i]);
                    for(int j=0; j<arr[i]; j++) {
                        clusters = product(clusters, getW(i+1));
                    }
                }
                Rational rFactor = new Rational(1, factor);
                for(int i=0; i<clusters.length; i++) {
                    clusters[i].setWeight(clusters[i].getWeight().times(rFactor));
                    list.add(clusters[i]);
                    ClusterDiagram fCluster = new ClusterDiagram(clusters[i]);
                    fCluster.addConnection(0, 1);
                    fCluster.addConnection(1, 0);
                    list.add(fCluster);
                }
            }
        }
        addEquivalents(list);
        h[n] = list.toArray(new ClusterDiagram[] {});
        return h[n];
    }
    
    public ClusterDiagram[] getC(int n) {
        if(n < c.length && c[n] != null) {
            return c[n];
        }
        if(n >= c.length) {
            c = Arrays.copyOf(c, n+1);
        }
        c[n] = difference(getH(n), getEta(n));
        return c[n];
    }

    public ClusterDiagram[] getB(int n) {
        if(n < b.length && b[n] != null) {
            return b[n];
        }
        if(approx == 0) {
            throw new RuntimeException("getB method not yet implemented without approximation");
        }
        return new ClusterDiagram[] {};
    }

    /**
     * Returns all diagrams of order n in the expansion for the indirect correlation function eta(r) = h(r) - c(r),
     * in the presence of any approximations currently in effect.  Evaluates eta via convolution of appropriate terms in
     * h(r) and c(r) up to order n-1.
     */
    public ClusterDiagram[] getEta(int n) {
        if(n < eta.length && eta[n] != null) {
            return eta[n];
        }
        if(n >= eta.length) {
            eta = Arrays.copyOf(eta, n+1);
        }
        LinkedList<ClusterDiagram> clusterList = new LinkedList<ClusterDiagram>();
        for(int i=0; i<n; i++) {
            int j = n - 1 - i;
            ClusterDiagram[] conv = convolution(getC(i), getH(j));
            for(int k=0; k<conv.length; k++) {
                clusterList.add(conv[k]);
            }
        }
        eta[n] = clusterList.toArray(new ClusterDiagram[] {});
        return eta[n];
    }
    
    /**
     * Returns all diagrams of order n in the sum [eta(r) + b(r)] (i.e., the sum of the indirect correlation
     * function and the bridge function). 
     */
    public ClusterDiagram[] getW(int n) {
        if(approx > 0) {
            return getEta(n);
        }
        if(n < w.length && w[n] != null) {
            return w[n];
        }
        if(n >= w.length) {
            w = Arrays.copyOf(w, n+1);
        }
        w[n] = sum(getEta(n), getB(n));
        return w[n];
    }
    
    /**
     * Returns an array of cluster diagrams with all e-bonds representing
     * the same sum of diagrams as the given f-bond cluster.
     */
    public static ArrayList<ClusterDiagram> toE(ClusterDiagram f) {
        ArrayList<ClusterDiagram> oldDiagrams = new ArrayList<ClusterDiagram>();
        oldDiagrams.add(f);
        ArrayList<ClusterDiagram> newDiagrams = new ArrayList<ClusterDiagram>();
        int totalConnections = f.getNumConnections();
        for (int iConnection=0; iConnection<totalConnections; iConnection++) {
            int node1 = -1;
            int node2 = -1;
            int connectionCount = 0;
            for (int i=0; i<f.mNumBody-1; i++) {
                for (int j=0; j<f.mNumBody; j++) {
                    if (f.mConnections[i][j] > i) {
                        if (connectionCount == iConnection) {
                            node1 = i;
                            node2 = f.mConnections[i][j];
                            break;
                        }
                        connectionCount++;
                    }
                    else if (f.mConnections[i][j] == -1) {
                        break;
                    }
                }
                if (node1 != -1 || node2 != -1) {
                    break;
                }
            }
            if (node1 == -1 || node2 == -1) {
                throw new RuntimeException("oops");
            }
            for (int i=0; i<oldDiagrams.size(); i++) {
                ClusterDiagram oldF = oldDiagrams.get(i);
                ClusterDiagram e0 = new ClusterDiagram(oldF);
                newDiagrams.add(new ClusterDiagram(e0));
                e0.deleteConnection(node1, node2);
                e0.deleteConnection(node2, node1);
                e0.setWeight(oldF.getWeight().times(new Rational(-1,1)));
                newDiagrams.add(e0);
            }
            
            oldDiagrams.clear();
            reduce(newDiagrams);
            ArrayList<ClusterDiagram> swapper = newDiagrams;
            newDiagrams = oldDiagrams;
            oldDiagrams = swapper;
        }
        return oldDiagrams;
    }
    
    /**
     * Returns an array of cluster diagrams with all e-bonds representing
     * the same sum of diagrams as the given f-bond clusters.
     */
    public static ClusterDiagram[] toE(ClusterDiagram[] f) {
        ArrayList<ClusterDiagram> e = new ArrayList<ClusterDiagram>();
        for (int i=0; i<f.length; i++) {
            e.addAll(toE(f[i]));
        }
        reduce(e);
        return e.toArray(new ClusterDiagram[e.size()]);
    }

    public static void main(String args[]) {
//        ClusterDiagram cluster1 = new ClusterDiagram(5, 2, Standard.ring(5));
//        cluster1.deleteConnection(0, 1);
//        cluster1.deleteConnection(1, 0);
//        ClusterDiagram cluster2 = new ClusterDiagram(4,2);
//        cluster2 = new ClusterDiagram(4, 2, Standard.ring(4));
//        cluster2.addConnection(0, 2);
//        cluster2.deleteConnection(0, 1);
//        cluster2.deleteConnection(1, 0);
//        ClusterDiagram product = ClusterOperations.convolution(cluster2, cluster1);
        
        //etomica.virial.junk.MyApplet applet = new etomica.virial.junk.MyApplet();
//        applet.init();
//        applet.starter.setDrawNumbersOfBlack(true);
//        applet.starter.setDrawNumbersOfWhite(true);
        
        int n = 3;
        ClusterOperations ops = new ClusterOperations();
        ClusterDiagram[] approxClusters = ops.getC(n);
//        applet.starter.addCluster(cluster1);
//        applet.starter.addCluster(cluster2);
//        applet.starter.addCluster(product);
        
        ArrayList<ClusterDiagram> list = new ArrayList<ClusterDiagram>();
        ClusterDiagram cluster = new ClusterDiagram(n+2, 2);
        ClusterGenerator gen = new ClusterGenerator(cluster);
        gen.setAllPermutations(false);
        gen.setOnlyConnected(false);
        gen.setOnlyDoublyConnected(true);
        gen.setExcludeArticulationPoint(true);
        gen.setExcludeArticulationPair(false);
        gen.setExcludeNodalPoint(true);
        gen.setMakeReeHover(false);
        //cluster.reset();
        gen.reset();
        cluster.setWeight(new Rational(1, cluster.mNumIdenticalPermutations));
        list.add(new ClusterDiagram(cluster));
        while(gen.advance()) {
            cluster.setWeight(new Rational(1, cluster.mNumIdenticalPermutations));
            list.add(new ClusterDiagram(cluster));
        }
        addEquivalents(list);
        ClusterDiagram[] trueClusters = list.toArray(new ClusterDiagram[] {});
        ClusterDiagram[] xs = difference(trueClusters, approxClusters);
        ClusterDiagram[] out = xs;
        out = integrate(out);
        out = integrate(out);
        out = makeReeHoover(out);
        for(int i=0; i<out.length; i++) {
            System.out.println(out[i]);
//            applet.starter.addCluster(out[i]);
        }
//        JPanel panel = applet.myPanel;
//        SimulationGraphic.makeAndDisplayFrame(panel, "ClusterOperation");
    }
    
    private static final ClusterDiagram zero = new ClusterDiagram(2, 2, new int[][] {}, new Rational(0,1));//cluster with zero weight
    private static final ClusterDiagram unity = new ClusterDiagram(2, 2, new int[][] {});//two root points with no bonds
    private static final ClusterDiagram f01 = new ClusterDiagram(2, 2);//two root points joined with a bond
    private static final ClusterDiagram f02f12 = new ClusterDiagram(3, 2, new int[][] {{0,2},{1,2}});
    private static final ClusterDiagram f01f02f12 = new ClusterDiagram(3, 2);
    private ClusterDiagram h[][] = new ClusterDiagram[][] {{f01}, {f02f12, f01f02f12}};
    private ClusterDiagram c[][] = new ClusterDiagram[][] {{f01}, {f01f02f12}};
    private ClusterDiagram eta[][] = new ClusterDiagram[][] {{zero}, {f02f12}};
    private ClusterDiagram w[][] = new ClusterDiagram[][] {eta[0], eta[1]};
    private ClusterDiagram b[][] = new ClusterDiagram[][] {{zero}, {zero}};
    public static final int NONE = 0;
    public static final int PY = 1;
    public static final int HNC = 2;
    private int approx = HNC;
}
