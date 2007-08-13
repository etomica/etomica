package etomica.virial.cluster;

import java.util.LinkedList;
import java.util.ListIterator;

import javax.swing.JPanel;

import etomica.graphics.SimulationGraphic;
import etomica.math.SpecialFunctions;
import etomica.math.discrete.FixedSumIterator;
import etomica.util.Arrays;
import etomica.util.Rational;

public class ClusterOperations {

    public ClusterOperations() {
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
        LinkedList list = new LinkedList();
        if(set1.length == 0 || set2.length == 0) {
            return new ClusterDiagram[] {};
        }
        for(int i=0; i<set1.length; i++) {
            for(int j=0; j<set2.length; j++) {
                list.add(product(set1[i],set2[j]));
            }
        }
        addEquivalents(list);
        return (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
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
        LinkedList list = new LinkedList();
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
        return (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
    }
    
    /**
     * Replaces all isomorphic clusters in the given list by a single cluster with weight
     * given by their sum.  Both the list and the the weights of the clusters it contains
     * may be altered by this method.
     */
    public static void addEquivalents(LinkedList list) {
        if(list.size() == 0 || list.size() == 1) {
            return;
        }
        ListIterator outer = list.listIterator(0);
        ClusterDiagram cluster1 = null;
        int j = 0;
        while(outer.hasNext()) {
            cluster1 = (ClusterDiagram)outer.next();
            ListIterator inner = list.listIterator(++j);
            while(inner.hasNext()) {
                ClusterDiagram cluster2 = (ClusterDiagram)inner.next();
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
                    } else {
                        cluster1.setWeight(cluster1.getWeight().plus(cluster2.getWeight()));
                        cluster2.setWeight(Rational.ZERO);
                    }
                }
            }
        }
        //remove clusters having zero weight
        ListIterator iterator = list.listIterator(0);
        while(iterator.hasNext()) {
            ClusterDiagram cluster = (ClusterDiagram)iterator.next();
            if(cluster.getWeight().numerator() == 0) {
                iterator.remove();
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
        LinkedList list = new LinkedList();
        for(int i=0; i<set1.length; i++) {
            list.add(new ClusterDiagram(set1[i]));
        }
        for(int i=0; i<set2.length; i++) {
            ClusterDiagram cluster2 = new ClusterDiagram(set2[i]);
            cluster2.setWeight(cluster2.getWeight().negate());
            list.add(cluster2);
        }
        addEquivalents(list);
        return (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
    }

    /**
     * Generates the set of diagrams obtained by adding the second set to the first set.
     * The returned array contains (copies of) all diagrams from the first set, with weights modified by adding
     * the weights of matching diagrams from the second set (if present), and all unmatched diagrams from
     * the second set. The returned array and all clusters in it are new instances.
     */
    public static ClusterDiagram[] sum(ClusterDiagram[] set1, ClusterDiagram[] set2) {
        LinkedList list = new LinkedList();
        for(int i=0; i<set1.length; i++) {
            list.add(new ClusterDiagram(set1[i]));
        }
        for(int i=0; i<set2.length; i++) {
            list.add(new ClusterDiagram(set2[i]));
        }
        addEquivalents(list);
        return (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
    }

    /**
     * Returns an array of clusters formed from taking the convolution of each of the diagrams
     * of the first set with each of the diagrams of the second set.  The given arrays and the
     * clusters they change are unchanged, and all clusters in the returned array are new.
     */
    public static ClusterDiagram[] convolution(ClusterDiagram[] set1, ClusterDiagram[] set2) {
        LinkedList list = new LinkedList();
        if(set1.length == 0 || set2.length == 0) {
            return new ClusterDiagram[] {};
        }
        for(int i=0; i<set1.length; i++) {
            for(int j=0; j<set2.length; j++) {
                list.add(convolution(set1[i],set2[j]));
            }
        }
        addEquivalents(list);
        return (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
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
                    rhClusters[i].addConnection(pairList[k][0], pairList[k][1]);
                }
            }
            if(neg) {//odd number of f-bonds introduced
                rhClusters[i].setWeight(rhClusters[i].getWeight().negate());
            }
        }
        return rhClusters;
    }
    
    public static ClusterDiagram[] makeReeHoover(ClusterDiagram[] clusters) {
        LinkedList list = new LinkedList();
        for(int i=0; i<clusters.length; i++) {
            ClusterDiagram[] rhClusters = makeReeHoover(clusters[i]);
            for(int j=0; j<rhClusters.length; j++) {
                list.add(rhClusters[j]);
            }
        }
        addEquivalents(list);
        return (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
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
            h = (ClusterDiagram[][])Arrays.resizeArray(h, n+1);
        }
        LinkedList list = new LinkedList();
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
        addEquivalents(list);
        h[n] =  (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
        return h[n];
    }
    
    public ClusterDiagram[] getC(int n) {
        if(n < c.length && c[n] != null) {
            return c[n];
        }
        if(n >= c.length) {
            c = (ClusterDiagram[][])Arrays.resizeArray(c, n+1);
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
     * in the presence of any approximations currently in effect.
     */
    public ClusterDiagram[] getEta(int n) {
        if(n < eta.length && eta[n] != null) {
            return eta[n];
        }
        if(n >= eta.length) {
            eta = (ClusterDiagram[][])Arrays.resizeArray(eta, n+1);
        }
        LinkedList clusterList = new LinkedList();
        for(int i=0; i<n; i++) {
            int j = n - 1 - i;
            ClusterDiagram[] conv = convolution(getC(i), getH(j));
            for(int k=0; k<conv.length; k++) {
                clusterList.add(conv[k]);
            }
        }
        eta[n] = (ClusterDiagram[])clusterList.toArray(new ClusterDiagram[] {});
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
            w = (ClusterDiagram[][])Arrays.resizeArray(w, n+1);
        }
        w[n] = sum(getEta(n), getB(n));
        return w[n];
    }
    
    public static void main(String args[]) {
////        ClusterDiagram cluster1 = new ClusterDiagram(5, 2, Standard.ring(5));
////        cluster1.deleteConnection(0, 1);
////        cluster1.deleteConnection(1, 0);
////        ClusterDiagram cluster2 = new ClusterDiagram(4,2);
////        cluster2 = new ClusterDiagram(4, 2, Standard.ring(4));
////        cluster2.addConnection(0, 2);
////        cluster2.deleteConnection(0, 1);
////        cluster2.deleteConnection(1, 0);
////        ClusterDiagram product = ClusterOperations.convolution(cluster2, cluster1);
//        
//        etomica.virial.junk.MyApplet applet = new etomica.virial.junk.MyApplet();
//        applet.init();
//        applet.starter.setDrawNumbersOfBlack(true);
//        applet.starter.setDrawNumbersOfWhite(true);
//        
//        int n = 4;
//        ClusterOperations ops = new ClusterOperations();
//        ClusterDiagram[] approxClusters = ops.getC(n);
////        applet.starter.addCluster(cluster1);
////        applet.starter.addCluster(cluster2);
////        applet.starter.addCluster(product);
//        
//        LinkedList list = new LinkedList();
//        ClusterDiagram cluster = new ClusterDiagram(n+2, 2);
//        ClusterGenerator gen = new ClusterGenerator(cluster);
//        gen.setAllPermutations(false);
//        gen.setOnlyConnected(false);
//        gen.setOnlyDoublyConnected(true);
//        gen.setExcludeArticulationPoint(true);
//        gen.setExcludeArticulationPair(false);
//        gen.setExcludeNodalPoint(true);
//        gen.setMakeReeHover(false);
//        //cluster.reset();
//        gen.reset();
//        cluster.setWeight(new Rational(1, cluster.mNumIdenticalPermutations));
//        list.add(new ClusterDiagram(cluster));
//        while(gen.advance()) {
//            cluster.setWeight(new Rational(1, cluster.mNumIdenticalPermutations));
//            list.add(new ClusterDiagram(cluster));
//        }
//        addEquivalents(list);
//        ClusterDiagram[] trueClusters = (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
//        ClusterDiagram[] xs = difference(trueClusters, approxClusters);
//        ClusterDiagram[] out = xs;
//        out = integrate(out);
//        out = integrate(out);
//        out = makeReeHoover(out);
//        for(int i=0; i<out.length; i++) {
//            applet.starter.addCluster(out[i]);
//        }
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
    private int approx = 1;//approximation level: 0 = none; 1 = PY; 2 = HNC
}
