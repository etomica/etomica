package etomica.virial.cluster;

import javax.swing.JPanel;

import etomica.graphics.SimulationGraphic;
//import etomica.virial.junk.MyApplet;

public class ClusterOperations {

    private ClusterOperations() {
    }
    
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
        for (int node1 = 0; node1 < cluster2.mNumBody; node1++) {
            int[] connections = cluster2.mConnections[node1];
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
    
    public static void main(String args[]) {
        ClusterDiagram cluster1 = new ClusterDiagram(5, 2, Standard.ring(5));
        cluster1.deleteConnection(0, 1);
        cluster1.deleteConnection(1, 0);
        ClusterDiagram cluster2 = new ClusterDiagram(5,2);
        cluster2 = new ClusterDiagram(5, 2, Standard.ring(5));
        ClusterDiagram product = ClusterOperations.product(cluster1, cluster2);
//        MyApplet applet = new MyApplet();
//        applet.init();
//        applet.starter.addCluster(cluster1);
//        applet.starter.addCluster(cluster2);
//        applet.starter.addCluster(product);
//        JPanel panel = applet.myPanel;
//        SimulationGraphic.makeAndDisplayFrame(panel, "ClusterOperation");
    }
}
