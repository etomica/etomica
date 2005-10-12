/*
 * Created on Sep 20, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial.cluster;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ClusterCollapsedTreeNode extends ClusterTreeNode {

    public ClusterCollapsedTreeNode(int numNodes, int numBondTypes) {
        super(numBondTypes);
        bonds = new int[numNodes][2];
    }
    
    public ClusterCollapsedTreeNode(ClusterTreeNode node1, ClusterTreeNode node2) {
        super(2);
        bonds = new int[node1.getNumNodes()+node2.getNumNodes()][2];
        int i;
        for (i=0; i<node1.getNumNodes(); i++) {
            bonds[i] = node1.getBond(i);
        }
        for (int j=i; j<bonds.length; j++) {
            bonds[j] = node2.getBond(j-i);
        }
        children[0] = node2.children[0];
        children[1] = node2.children[1];
        coefficient = node2.coefficient;
    }
    
    public int getNumNodes() {
        return bonds.length;
    }
    
    public int[] getBond(int i) {
        return bonds[i];
    }
    
    public void setBonds(int [][]newBonds) {
        for (int i=0; i<bonds.length; i++) {
            bonds[i][0] = newBonds[i][0];
            bonds[i][1] = newBonds[i][1];
        }
    }
    
    public double value(double[][] fValues) {
        double p = 1;
        for (int i=0; i<bonds.length; i++) {
            int[] iBond = bonds[i];
            int bondNum = iBond[0];
            int fNum = iBond[1];
            p *= fValues[bondNum][fNum];
            if (Double.isNaN(p)) {
                System.out.println("oops8");
            }
        }
        if (Double.isNaN(p)) {
            System.out.println("oops8");
        }
        return p;
    }
    
    private final int[][] bonds;
}
