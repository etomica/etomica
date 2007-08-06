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
     * the clusters it contains are unchanged.
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
     * the clusters in them are not changed.
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
     * Replaces all isomorphic clusters in the given list by a single cluster with weight
     * given by their sum.  Both the list and the the weights of the clusters it contains
     * may be altered by this method.
     */
    public static void addEquivalents(LinkedList list) {
        ListIterator outer = list.listIterator(0);
        int j = 0;
        while(outer.hasNext()) {
            ClusterDiagram cluster1 = (ClusterDiagram)outer.next();
            ListIterator inner = list.listIterator(++j);
            while(inner.hasNext()) {
                ClusterDiagram cluster2 = (ClusterDiagram)inner.next();
                if(cluster1.equals(cluster2)) {
                    cluster1.setWeight(cluster1.getWeight().plus(cluster2.getWeight()));
                    cluster2.setWeight(Rational.ZERO);
                }
            }
        }
        //removed clusters having zero weight
        while(outer.hasPrevious()) {
            ClusterDiagram cluster = (ClusterDiagram)outer.previous();
            if(cluster.getWeight().numerator() == 0) {
                outer.remove();
            }
        }
    }

    
    /**
     * Returns a cluster that is a convolution of the given cluster, which is obtained by joining a root point
     * from each diagram and turning it into a field point.  Each given cluster must have exactly 2 root points.
     * The weight of the returned cluster is the product of the weights of the given clusters.
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
     * the seconds set with their weights negated.
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
     * Returns an array of clusters formed from taking the convolution of each of the diagrams
     * of the first set with each of the diagrams of the second set.
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
                list.add(fCluster);
            }
        }
        addEquivalents(list);
        return (ClusterDiagram[])list.toArray(new ClusterDiagram[] {});
    }
    
    public ClusterDiagram[] getC(int n) {
        if(n < c.length && c[n] != null) {
            return c[n];
        }
        return null;
    }

    public ClusterDiagram[] getB(int n) {
        if(n < b.length && b[n] != null) {
            return b[n];
        }
        return null;
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
        return (ClusterDiagram[])clusterList.toArray(new ClusterDiagram[] {});
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
        return null;
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
//        ClusterOperations ops = new ClusterOperations();
//        ClusterDiagram[] clusters = ops.getH(2);
////        applet.starter.addCluster(cluster1);
////        applet.starter.addCluster(cluster2);
////        applet.starter.addCluster(product);
//        
//        for(int i=0; i<clusters.length; i++) {
//            applet.starter.addCluster(clusters[i]);
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
