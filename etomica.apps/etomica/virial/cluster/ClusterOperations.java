package etomica.virial.cluster;

import javax.swing.JPanel;

import etomica.graphics.SimulationGraphic;

public class ClusterOperations {

    private ClusterOperations() {
    }
    
    /**
     * Returns a new cluster that is the product of the given clusters.  Input clusters must have the
     * same number of root points and my not both have bonds between the same root points.
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
        return productCluster; 
    }
    
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
        return convolutionCluster;
    }
    
    public static void main(String args[]) {
        ClusterDiagram cluster1 = new ClusterDiagram(5, 2, Standard.ring(5));
        cluster1.deleteConnection(0, 1);
        cluster1.deleteConnection(1, 0);
        ClusterDiagram cluster2 = new ClusterDiagram(4,2);
//        cluster2 = new ClusterDiagram(4, 2, Standard.ring(4));
//        cluster2.addConnection(0, 2);
//        cluster2.deleteConnection(0, 1);
//        cluster2.deleteConnection(1, 0);
        ClusterDiagram product = ClusterOperations.convolution(cluster2, cluster1);
//        etomica.virial.junk.MyApplet applet = new etomica.virial.junk.MyApplet();
//        applet.init();
//        applet.starter.setDrawNumbersOfBlack(true);
//        applet.starter.setDrawNumbersOfWhite(true);
//        applet.starter.addCluster(cluster1);
//        applet.starter.addCluster(cluster2);
//        applet.starter.addCluster(product);
//        JPanel panel = applet.myPanel;
//        SimulationGraphic.makeAndDisplayFrame(panel, "ClusterOperation");
    }
}
