/*
 * Created on Sep 19, 2005
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
public class ClusterTreeNode {

    public ClusterTreeNode(int n) {
        children = new ClusterTreeNode[n];
    }
    
    public int getNumNodes() {
        return 1;
    }
    
    public int[] getBond(int i) {
        return new int[]{bondNum,fNum};
    }
    
    public void setBonds(int[] bonds) {
        bondNum = bonds[0];
        fNum = bonds[1];
    }
    
    public double value(double[][] fValues) {
        return fValues[bondNum][fNum];
    }
    
    public final ClusterTreeNode[] children;
    public double coefficient;
    private int bondNum;
    private int fNum;
}
