/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.util.random.IRandom;

/**
 * Holds information about a cluster diagram, including the bonds, root points,
 * and score.
 * @author andrew
 */
public class ClusterDiagramTree implements java.io.Serializable, DiagramSource {

    /**
     * Constructs a cluster having numBody points and numRootPoints root points.
     */
    public ClusterDiagramTree(int numBody, int numBondTypes) {
        mNumBody = numBody;
        mNumBondTypes = numBondTypes;
        mRootNode = new ClusterTreeNode(numBondTypes);
        fValue = new double[mNumBody*(mNumBody+1)/2];
        fValue[0] = 1.0;
        iBond = new int[mNumBody*(mNumBody+1)/2];
        currentNode = new ClusterTreeNode[numBody*(mNumBody-1)/2];
    }
    
    public int pointCount() {
        return mNumBody;
    }
    
    public void setPrefactor(double f) {
        prefactor = f;
    }
    
    public void setHistogramming(boolean flag) {
        histogramming = flag;
        if (histogramming) {
            values = new double[numDiagrams];
        }
    }
    
    public void setDiagramSampling(boolean flag) {
        doDiagramSampling = flag;
        if (doDiagramSampling) {
            values = new double[numDiagrams];
        }
    }
    
    public void setRandom(IRandom newRandom) {
        random = newRandom;
    }
    
    /**
     * Adds the given diagram to the tree, adding nodes as needed
     */
    public void addClusterDiagram(ClusterDiagram newDiagram) {
        ClusterTreeNode node = mRootNode;
        int bondCount = 0;
        for (int i=0; i<mNumBody-1; i++) {
            int k=0;
            for (int j=i+1; j<mNumBody; j++) {
                int l = newDiagram.mConnections[i][k];
                while (l < j && l != -1) {
                    // skip bonds pointing back at previous points
                    l = newDiagram.mConnections[i][++k];
                }
                if (l == -1) {
                    // no more bonds
                    if (mNumBondTypes == 2) {
                        // add any remaining points as being of the second bond type 
                        if (node.children[1] == null) {
                            node.children[1] = new ClusterTreeNode(mNumBondTypes);
                            node.children[1].setBonds(new int[]{bondCount,1});
                        }
                        bondCount++;
                        node = node.children[1];
                        continue;
                    }
                    // go to the next i point
                    break;
                }
                if (mNumBondTypes == 2 && l > j) {
                    // add any points not lists in mConnections as the second bond type 
                    if (node.children[1] == null) {
                        node.children[1] = new ClusterTreeNode(mNumBondTypes);
                        node.children[1].setBonds(new int[]{bondCount,1});
                    }
                    bondCount++;
                    node = node.children[1];
                }
                    
                if (l == j) {
                    // bond between i and j exists
                    if (node.children[0] == null) {
                        node.children[0] = new ClusterTreeNode(mNumBondTypes);
                        node.children[0].setBonds(new int[]{bondCount,0});
                    }
                    bondCount++;
                    node = node.children[0];
                    k++;
                }
            }
        }
        node.coefficient = newDiagram.mReeHooverFactor;
        numDiagrams++;
    }
    
    public int getNDiagrams() {
        return numDiagrams;
    }
    
    public double[] getLastValues() {
        return values;
    }

    public void writeDiagrams() {
        ClusterTreeNode node = mRootNode;
        String[] bonds = new String[mNumBody*(mNumBody+1)/2];
        bonds[0] = "";
        for (int i=0; i<iBond.length; i++) {
            iBond[i] = -1;
        }
        int nBonds = 0; //which bond we're on now (0,1) vs. (0,2) etc.  nBonds=0 means (0,1)
        int point1 = 0;
        int point2 = 1;
        while (true) {
            if (node.coefficient != 0) {
                System.out.println(node.coefficient+"   "+bonds[nBonds]);
            }
outer:      while (true) {
                while (iBond[nBonds] < mNumBondTypes-1) {
                    iBond[nBonds]++;
                    if (node.children[iBond[nBonds]] != null) {
                        node = node.children[iBond[nBonds]];
                        break outer;
                    }
                }
                // we made it through all of the bonds for this node.
                // now go back up.
//                node = node.parent;
                if (node == null) {
                    // back up to the root node and it has no more children
                    return;
                }
                // reset bond type counter for this level
                iBond[nBonds] = -1;
                nBonds--;
                point2--;
                if (point2 == point1) {
                    point1--;
                    point2 = mNumBody-1;
                }
            }
            nBonds++;
            bonds[nBonds] = bonds[nBonds-1]+"   "+iBond[nBonds-1]+" "+point1+" "+point2;
            point2++;
            if (point2 == mNumBody) {
                point1++;
                point2 = point1+1;
            }
        }
    }
    
    public void collapse() {
        ClusterTreeNode node = mRootNode;
        for (int i=0; i<iBond.length; i++) {
            iBond[i] = -1;
        }
        int nBonds = 0; //which bond we're on now (0,1) vs. (0,2) etc.  nBonds=0 means (0,1)
        while (true) {
            while (true) {
                //assume mNumBondTypes == 2
                ClusterTreeNode node2 = null;
                if (node.children[0] == null) {
                    if (node.children[1] == null) {
                        // we made it through all of the bonds for this node.
                        // now go back up.
                        if (nBonds == 0) {
                            // back up to the root node and it has no more children
                            return;
                        }
                        // reset bond type counter for this level
                        iBond[nBonds] = -1;
                        nBonds--;
                        node = currentNode[nBonds];
                        continue;
                    }
                    node2 = node.children[1];
                }
                else if (node.children[1] == null) {
                    node2 = node.children[0];
                }
                if (node2 != null) {
                    ClusterCollapsedTreeNode mergedNode = new ClusterCollapsedTreeNode(node,node2);
                    ClusterTreeNode parent = currentNode[nBonds-1];
                    if (parent.children[0] == node) {
                        parent.children[0] = mergedNode;
                    }
                    else {
                        parent.children[1] = mergedNode;
                    }
                    node = mergedNode;
                    continue;
                }
                iBond[nBonds]++;
                if (iBond[nBonds] == mNumBondTypes) {
                    // we made it through all of the bonds for this node.
                    // now go back up.
                    if (nBonds == 0) {
                        // back up to the root node and it has no more children
                        return;
                    }
                    // reset bond type counter for this level
                    iBond[nBonds] = -1;
                    nBonds--;
                    node = currentNode[nBonds];
                    continue;
                }
                break;
            }
            currentNode[nBonds] = node;
            node = node.children[iBond[nBonds]];
            nBonds++;
        }
    }
    
    public double value(double[][] fValues) {
        double sum = 0.0;
        int iDiagram = 0;
        if (histogramming) {
            for (int i=0; i<values.length; i++) {
                values[i] = 0;
            }
        }
        ClusterTreeNode node = mRootNode;
        for (int i=0; i<iBond.length; i++) {
            iBond[i] = -1;
        }
        int nBonds = 0; //which bond we're on now (0,1) vs. (0,2) etc.  nBonds=0 means (0,1)
        double sumAbs = 0;
        while (true) {
            if (node.coefficient != 0) {
                if (histogramming || doDiagramSampling) {
                    values[iDiagram++] = node.coefficient * fValue[nBonds];
                }
                if (doDiagramSampling) {
//                    System.out.print(node.coefficient * fValue[nBonds]+" ");
                    sumAbs += Math.abs(node.coefficient * fValue[nBonds]);
                }
                sum += node.coefficient * fValue[nBonds];
                if (nBonds == 0) {
                    // can only get in here if root node is the only node.
                    return prefactor * sum;
                }
                // reset bond type counter for this level
                iBond[nBonds] = -1;
                nBonds--;
                // we made it through all of the bonds for this node.
                // now go back up.
                node = currentNode[nBonds];
            }
            double v;
outer:      while (true) {
                while (iBond[nBonds] < mNumBondTypes-1) {
                    iBond[nBonds]++;
                    ClusterTreeNode childNode = node.children[iBond[nBonds]];
                    v = childNode.value(fValues);
                    if (v != 0) {
                        currentNode[nBonds] = node;
                        node = childNode;
                        break outer;
                    }
                }
                // we made it through all of the bonds for this node.
                // now go back up.
                if (nBonds == 0) {
                    if (doDiagramSampling) {
                        double r = random.nextDouble()*sumAbs;
                        double sum2 = 0;
                        for (pickedDiagram = 0; pickedDiagram<numDiagrams; pickedDiagram++) {
                            sum2 += Math.abs(values[pickedDiagram]);
                            if (sum2 > r) {
//                                System.out.println("==> "+values[pickedDiagram]);
                                return numDiagrams*values[pickedDiagram];
                            }
                        }
//                        System.out.println("=oops=> "+values[pickedDiagram-1]);
                        pickedDiagram = numDiagrams-1; 
                        return numDiagrams*values[pickedDiagram];
                    }
//                    System.exit(0);
                    // back up to the root node and it has no more children
                    if (Double.isNaN(prefactor*sum)) {
                        System.out.println("oops4");
                    }
                    return prefactor*sum;
                }
                // reset bond type counter for this level
                iBond[nBonds] = -1;
                nBonds--;
                node = currentNode[nBonds];
            }
            fValue[nBonds+1] = fValue[nBonds] * v;
            nBonds++;
        }
    }
    
    public int getPickedDiagram() {
        return pickedDiagram;
    }

    private static final long serialVersionUID = 1L;
    protected IRandom random;
    private final int mNumBody;
    private final int mNumBondTypes;
    private final ClusterTreeNode mRootNode;
    private final double[] fValue;
    private double prefactor;
    private int[] iBond;
    private int numDiagrams = 0;
    private double[] values;
    private boolean histogramming = false;
    private boolean doDiagramSampling = false;
    private int pickedDiagram = -1;
    private final ClusterTreeNode[] currentNode;
}
