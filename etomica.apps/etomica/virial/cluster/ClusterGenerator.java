package etomica.virial.cluster;

import java.io.FileWriter;
import java.io.IOException;


/**
 * Generates all clusters matching certain criteria.
 * @author andrew
 */
public class ClusterGenerator implements java.io.Serializable {
    
    public ClusterGenerator(ClusterDiagram aCluster) {
        mCluster = aCluster;
        excludeRootPermutations = true;
        mStack = new int[aCluster.mNumBody];
        mMarked = new boolean[aCluster.mNumBody];
        // use an upper bound for # of connections
        int nConnections = mCluster.mNumBody * (mCluster.mNumBody-1);
        mConnectList = new int[nConnections][2];
        mConnectStack = new int[nConnections];
        mStartConnection = 0;
        mStopConnection = Integer.MAX_VALUE;
    }

    /**
     * Set whether to only generate clusters that are singly-connected. 
     */
    public void setOnlyConnected(boolean flag) {
        mOnlyConnected = flag;
    }

    /**
     * Set whether to only generate clusters that are doubly-connected.
     */
    public void setOnlyDoublyConnected(boolean flag) {
        mOnlyDoublyConnected = flag;
    }

    /**
     * Set whether to exclude clusters having a nodal point.
     */
    public void setExcludeNodalPoint(boolean flag) {
        mExcludeNodalPoint = flag;
    }

    /**
     * Set whether to exclude clusters having an articulation point.
     */
    public void setExcludeArticulationPoint(boolean flag) {
        mExcludeArticulationPoint = flag;
    }

    /**
     * Set whether to exclude clusters having an articulation pair.
     */
    public void setExcludeArticulationPair(boolean flag) {
        mExcludeArticulationPair = flag;
    }

    /**
     * Set whether to generate all permutations.  If false, only
     * one cluster from each set of permutations will be generated.  
     * That cluster has the greatest "score" among its set of permutations,
     * and the number of identical permutations will be calculated and stored
     * in the mNumIdenticalPermutations field. 
     */
    public void setAllPermutations(boolean flag) {
        mAllPermutations = flag;
    }

    /**
     * Set whether to generate Ree-Hoover diagrams.  If set, then
     * connected points have an f-bond and all others have an e-bond.
     */
    public void setMakeReeHover(boolean flag) {
        mReeHoover = flag;
        if (flag) {
            mClusterCopy = new ClusterDiagram(mCluster);
            mReeHooverGenerator = new ClusterGenerator(mClusterCopy);
            mReeHooverGenerator.setAllPermutations(true);
            mReeHooverGenerator.setOnlyDoublyConnected(true);
        }
    }

    public static boolean swapAllAndCompare(boolean excludeRootPermutations, final ClusterDiagram cluster, int[] score) {
        return swapAllAndCompare(0, cluster.mNumBody-1, excludeRootPermutations, cluster, score);
    }

    /**
     * Determines whether any cluster permuted from the given one by
     * exchanging points from startPoint to stopPoint (inclusive) has a higher
     * score than that passed in. If such a cluster is located, the method returns
     * true immediately. The method always restores the connections to their
     * state when the method was called.  If the permutations of the given cluster 
     * do not have a higher score, the cluster's mNumIdenticalPermutations field is 
     * updated to include the number of identical permutations encountered while 
     * swapping.   
     */
    private static boolean swapAllAndCompare(int startPoint, int stopPoint, boolean excludeRootPermutations,
            final ClusterDiagram cluster, int[] score) {
        // recursively swap everything between startPoint+1 and stopPoint
        // this is basically the 0th iteration below
        if (startPoint+1 < stopPoint && swapAllAndCompare(startPoint + 1, stopPoint, excludeRootPermutations, cluster, score)) {
            return true;
        }
        if(cluster.isRootPoint(startPoint) && excludeRootPermutations) {
            return false;
        }
        // recursively swap startPoint with startPoint+1 to stopPoint
        for (int i = startPoint + 1; i < stopPoint + 1; i++) {
            if (cluster.isRootPoint(startPoint) == cluster.isRootPoint(i)) {
                cluster.swap(startPoint, i);
                boolean isGreaterThan = cluster.scoreGreaterThan(score);
                isGreaterThan = isGreaterThan || swapAllAndCompare(startPoint + 1, stopPoint, excludeRootPermutations, cluster, score);
                cluster.swap(startPoint, i);
                if (isGreaterThan) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Determines whether any cluster permuted from the given one by
     * exchanging points from startPoint to stopPoint (inclusive) has a higher
     * score than that passed in. If such a cluster is located, the method returns
     * true immediately. The method always restores the connections to their
     * state when the method was called.  If the permutations of the given cluster 
     * do not have a higher score, the cluster's mNumIdenticalPermutations field is 
     * updated to include the number of identical permutations encountered while 
     * swapping.   
     */
    public static void findMaxScore(boolean excludeRootPermutations, final ClusterDiagram cluster, int[] score) {
        cluster.calcScore(score);
        findMaxScore(0, cluster.mNumBody-1, excludeRootPermutations, cluster, score);
    }
    
    private static void findMaxScore(int startPoint, int stopPoint, boolean excludeRootPermutations, 
            final ClusterDiagram cluster, int[] score) {
        // recursively swap everything between startPoint+1 and stopPoint
        // this is basically the 0th iteration below
        if (startPoint+1 < stopPoint) {
            findMaxScore(startPoint + 1, stopPoint, excludeRootPermutations, cluster, score);
        }
        if(cluster.isRootPoint(startPoint) && excludeRootPermutations) {
            return;
        }
        // recursively swap startPoint with startPoint+1 to stopPoint
        for (int i = startPoint + 1; i < stopPoint + 1; i++) {
            if (cluster.isRootPoint(startPoint) == cluster.isRootPoint(i)) {
                cluster.swap(startPoint, i);
                if (cluster.scoreGreaterThan(score)) {
                    cluster.calcScore(score);
                }
                findMaxScore(startPoint + 1, stopPoint, excludeRootPermutations, cluster, score);
                cluster.swap(startPoint, i);
            }
        }
    }

    /**
     * Determines if the cluster is connected after removal of the
     * given node.
     */
    public boolean isConnected(int exclude) {
        for (int i = 0; i < mCluster.mNumBody; i++) {
            mMarked[i] = false;
        }
        mMarked[exclude] = true;
        int numMarked = 1;
        mStack[0] = 0;
        if (exclude == 0) {
            mStack[0] = 1;
            numMarked = 2;
            mMarked[1] = true;
        }
        mStack[1] = -1;
        return isConnectedReally(numMarked);
    }

    /**
     * Determines if the cluster is connected.
     */
    public boolean isConnected() {
        for (int i = 0; i < mCluster.mNumBody; i++) {
            mMarked[i] = false;
        }
        mMarked[0] = true;
        mStack[0] = 0;
        mStack[1] = -1;
        return isConnectedReally(1);
    }
    
    /**
     * Does the real work for the public is*Connected methods.  The
     * method walks through the cluster ensuring that every point is 
     * reachable.
     */
    private boolean isConnectedReally(int numMarked) {
        int stackLength = 1;
        while (stackLength > 0) {
            // pop the mStack
            int thisNode = mStack[stackLength - 1];
            stackLength--;
            for (int i = 0;; i++) {
                int connectedNode = mCluster.mConnections[thisNode][i];
                if (connectedNode == -1) {
                    // end of connections for thisNode
                    break;
                }
                if (!mMarked[connectedNode]) {
                    mMarked[connectedNode] = true;
                    numMarked++;
                    if (numMarked == mCluster.mNumBody) {
                        return true;
                    }
                    mStack[stackLength] = connectedNode;
                    stackLength++;
                }
            }
        }
        return false;
    }

    /**
     * Convenience method to determine whether the cluster is elementary or not.
     */
    public boolean isElementary() {
        for (int i = mCluster.mNumBody - 1; i >= 0; i--) {
            if (isNodalPoint(i) || isArticulationPoint(i)) return false;
            for (int j = i - 1; j >= 0; j--) {
                if (isArticulationPair(i, j)) return false;
            }
        }
        return true;
    }

    /**
     * Determines if the given node is a nodal point.
     */
    public boolean isNodalPoint(int exclude) {
        if (mCluster.isRootPoint(exclude) || mCluster.getNumRootPoints() < 2) {
            return false;
        }
        int stackLength = 1;
        int i;
        int thisNode, connectedNode;
        int rootPointsFound = 1;

        for (i = 0; i < mCluster.mNumBody; i++) {
            mMarked[i] = false;
        }
        mMarked[exclude] = true;
        mStack[0] = 0;
        mStack[1] = -1; // -1 terminates the mStack
        mMarked[0] = true;

        while (stackLength > 0) {
            // pop the mStack
            thisNode = mStack[stackLength - 1];
            stackLength--;
            for (i = 0;; i++) {
                connectedNode = mCluster.mConnections[thisNode][i];
                if (connectedNode == -1) {
                    // end of connections for thisNode
                    break;
                }
                if (!mMarked[connectedNode]) {
                    if (mCluster.isRootPoint(connectedNode)) {
                        if (++rootPointsFound == mCluster.getNumRootPoints()) {
                            return false;
                        }
                    }
                    mMarked[connectedNode] = true;
                    mStack[stackLength] = connectedNode;
                    stackLength++;
                }
            }
        }
        return true;
    }

    /**
     * Determines if the given point is an articulation point.
     */
    public boolean isArticulationPoint(int exclude) {
        int i;
        int numMarksNeeded;
        // # of root points that haven't been excluded
        int nRootPoints = mCluster.getNumRootPoints();

        for (i = 0; i < mCluster.mNumBody; i++) {
            mMarked[i] = false;
        }
        mMarked[exclude] = true;
        numMarksNeeded = mCluster.mNumBody - mCluster.getNumRootPoints();
        if (!mCluster.isRootPoint(exclude)) {
            numMarksNeeded--;
        }
        else {
            // excluded a root point, so don't make isArticulated try to find it
            nRootPoints--;
        }
        mStack[0] = mCluster.getNumRootPoints();
        while (mMarked[mStack[0]]) {
            mStack[0]++;
        }
        mMarked[mStack[0]] = true;
        mStack[1] = -1; // -1 terminates the mStack

        return isArticulated(numMarksNeeded,nRootPoints);
    }

    /**
     * determines if the given pair of points is an articulation pair.
     */
    public boolean isArticulationPair(int exclude1, int exclude2) {
        int i;
        int numMarksNeeded;
        // # of root points that haven't been excluded
        int nRootPoints = mCluster.getNumRootPoints();

        for (i = 0; i < mCluster.mNumBody; i++) {
            mMarked[i] = false;
        }
        // numMarked only includes marked filled points
        numMarksNeeded = mCluster.mNumBody - mCluster.getNumRootPoints();
        mMarked[exclude1] = true;
        if (!mCluster.isRootPoint(exclude1)) {
            numMarksNeeded--;
        }
        else {
            // excluded a root point, so don't make isArticulated try to find it
            nRootPoints--;
        }
        mMarked[exclude2] = true;
        if (!mCluster.isRootPoint(exclude2)) {
            numMarksNeeded--;
        }
        else {
            nRootPoints--;
        }
        mStack[0] = mCluster.getNumRootPoints();
        while (mMarked[mStack[0]]) {
            mStack[0]++;
            if (mStack[0] == mCluster.mNumBody) {
                // excluded all filled points
                return false;
            }
        }
        mMarked[mStack[0]] = true;
        mStack[1] = -1; // -1 terminates the mStack
        
        return isArticulated(numMarksNeeded, nRootPoints);
    }

    /**
     * Method that does the real work for the isArticulation* methods.
     * The method walks through the cluster looking for root nodes.  If
     * it finds a root point in every section, it returns true.  If all
     * root points were deleted in calling methods, the it returns true if
     * all filled points are in one section.
     */
    private boolean isArticulated(int numMarksNeeded, int nRootPoints) {
        int stackLength = 1;
        int numMarked = 1;
        boolean reachedReferenceNode;
        int pass = 1;

        while (true) {
            reachedReferenceNode = false;
            while (stackLength > 0) {
                // pop the mStack
                int thisNode = mStack[stackLength - 1];
                stackLength--;
                for (int i = 0;; i++) {
                    int connectedNode = mCluster.mConnections[thisNode][i];
                    if (connectedNode == -1) {
                        // end of connections for thisNode
                        break;
                    }
                    if (!mMarked[connectedNode]) {
                        if (!mCluster.isRootPoint(connectedNode)) {
                            numMarked++;
                            if (numMarked == numMarksNeeded
                                    && reachedReferenceNode) {
                                return false;
                            }
                        }
                        else {
                            if (numMarked == numMarksNeeded) {
                                // we've marked all the filled points and found
                                // a reference node
                                return false;
                            }
                            reachedReferenceNode = true;
                        }
                        mMarked[connectedNode] = true;
                        mStack[stackLength] = connectedNode;
                        stackLength++;
                    }
                }
            }
            if (numMarksNeeded == numMarked && nRootPoints == 0) {
                // all reference nodes were deleted and all filled points were
                // marked
                return false;
            }
            // If we didn't reach a reference node, the cluster must be articulated since
            // we didn't delete all reference nodes.
            // If this is the nth pass (n=# of root points) but we haven't marked everything
            // then there must be another chunk out there that doesn't have a reference node
            if (!reachedReferenceNode || pass == nRootPoints) {
                return true;
            }
            pass++;
            // find an unmarked filled point
            int thisNode = mCluster.getNumRootPoints()+1;
            while (mMarked[thisNode]) {
                thisNode++;
            }
            mMarked[thisNode] = true;
            numMarked++;
            mStack[0] = thisNode;
            stackLength = 1;
        }
    }

    /**
     * Determines if the cluster has a maximum score. If the score
     * can be increased by swapping filled points or by swapping reference
     * nodes, then one of its permutations must have a higher score and the
     * method returns false.  If the cluster does have a maximum score, the 
     * mNumIdenticalPermutations is calculated as 1 + the number of permuted 
     * clusters with the same score as the unpermuted cluster.
     */
    public boolean isMaximumScore() {
        int[] score = new int[mCluster.mNumBody/2+1];
        mCluster.calcScore(score);
        mCluster.mNumIdenticalPermutations = 1;
        return !swapAllAndCompare(excludeRootPermutations, mCluster, score);
    }

    /**
     * Calculates the Ree-Hoover coefficient (coefficient for Ree-Hoover 
     * diagram / coefficient for Mayer diagram) for mCluster and sets mCluster's
     * mReeHooverFactor to it.
     */
    public void calcReeHoover() {
        if (!mOnlyDoublyConnected) {
            throw new IllegalStateException("Ree-Hoover diagrams only work for doubly-connected cluster diagrams");
        }
        int originalConnections = mCluster.getNumConnections();
        mCluster.copyTo(mClusterCopy);
        mReeHooverGenerator.reset();
        mCluster.mReeHooverFactor = 1;
        while (mReeHooverGenerator.advance()) {
            int numRemoved = originalConnections - mClusterCopy.getNumConnections();
            mCluster.mReeHooverFactor += (numRemoved % 2 == 0) ? 1 : -1;
        }
    }
    
    /**
     * Resets the cluster generator.  This resets the cluster, and prepares the
     * generator to generate clusters matching the appropriate criteria.
     */
    public void reset() {
        for (int k = 0, i = 0; i < mCluster.mNumBody; i++) {
            int[] iConnections = mCluster.mConnections[i];
            for (int j = 0; iConnections[j] != -1; j++) {
                if (i < iConnections[j]) {
                    mConnectList[k][0] = i;
                    mConnectList[k][1] = iConnections[j];
                    k++;
                }
            }
        }
        
        // get # of connections in fully-connected cluster (should be k)
        mTotalNumConnections = mCluster.getNumConnections();
        mCurrentConnection = mStartConnection;
        mConnectStackLen = 0;
        if (!mAllPermutations && !isMaximumScore()) {
            throw new RuntimeException("initial cluster not maxed");
        }
    }

    public void setStartConnection(int i) {
        mStartConnection = i;
    }
    
    public void setStopConnection(int i) {
        mStopConnection = i;
    }
    
    /**
     * Advances the cluster to the next cluster matching the given criteria.  Returns
     * false if no more clusters exist. 
     */
    public boolean advance() {
        for (; mCurrentConnection < mTotalNumConnections || mConnectStackLen > 0; mCurrentConnection++) {
            if (mCurrentConnection == mTotalNumConnections) {
                // we made it through the whole list. pop a connection off the
                // mStack and go through the list again from there.
                mCurrentConnection = mConnectStack[mConnectStackLen - 1];
                mConnectStackLen--;
                int node1 = mConnectList[mCurrentConnection][0];
                int node2 = mConnectList[mCurrentConnection][1];
//                System.out.println("adding back connection between "+node1+" and "+node2);
                mCluster.addConnection(node1, node2);
                mCluster.addConnection(node2, node1);
                continue;
            }
            int node1 = mConnectList[mCurrentConnection][0];
            int node2 = mConnectList[mCurrentConnection][1];
            if (mConnectStackLen < 1) {
//              System.out.println(mConnectStackLen+" "+node1+" "+node2);
                if (mCurrentConnection == mStopConnection)
                    return false;
            }
//            System.out.println("deleting connection "+mCurrentConnection+" betweeen "+node1+" and "+node2);
            if (!mCluster.deleteConnection(node1, node2)) {
                throw new RuntimeException(
                        "failed to delete a connection that should have existed");
            }
            if (!mCluster.deleteConnection(node2, node1)) {
                throw new RuntimeException(
                        "failed to delete a connection that should have existed");
            }
            
            if (meetsAllCriteria()) {
                mConnectStack[mConnectStackLen] = mCurrentConnection;
                mConnectStackLen++;
                if (mAllPermutations || isMaximumScore()) {
                    if (mReeHoover) {
                        calcReeHoover();
                        if (mCluster.mReeHooverFactor == 0) {
                            // cluster diagram is cancelled out in ReeHoover formulation
                            continue;
                        }
                    }
                    mCluster.sort();
                    mCurrentConnection++;
                    return true;
                }
            }
            else {
                mCluster.addConnection(node1, node2);
                mCluster.addConnection(node2, node1);
            }
        }
        return false;
    }

    //returns true if the current cluster meets all indicated topological criteria
    private boolean meetsAllCriteria() {
        if (mExcludeNodalPoint) {
            for (int i = mCluster.mNumBody - 1; i >= 0; i--) {
                if (isNodalPoint(i)) {
                    return false;
                }
            }
        }
        if (mExcludeArticulationPoint) {
            for (int i = mCluster.mNumBody - 1; i >= 0; i--) {
                if (isArticulationPoint(i)) {
                    return false;
                }
            }
        }
        if (mExcludeArticulationPair) {
            for (int i = mCluster.mNumBody - 1; i >= 0; i--) {
                for (int j = i - 1; j >= 0; j--) {
                    if (isArticulationPair(i, j)) {
                        return false;
                    }
                }
            }
        }
        if (mOnlyDoublyConnected) {
            for (int i = mCluster.mNumBody - 1; i >= 0; i--) {
                if (!isConnected(i)) {
                    return false;
                }
            }
        }
        if (mOnlyConnected) {
            if(!isConnected()) {
                return false;
            }
        }
        return true;
    }

    /**
     * Demonstrates how to use this class.
     */
    public static void main(String[] args) {
        int nBody = 9;
        int nRootPoints = 0;
        if (args.length > 0) nBody = Integer.parseInt(args[0]);
        if (args.length > 1) nRootPoints = Integer.parseInt(args[1]);
        ClusterDiagram cluster = new ClusterDiagram(nBody,nRootPoints);
        ClusterGenerator generator = new ClusterGenerator(cluster);
        generator.setAllPermutations(false);
        generator.setOnlyDoublyConnected(nRootPoints == 0);
        generator.setExcludeArticulationPoint(nRootPoints > 0);
        generator.setExcludeArticulationPair(nRootPoints > 0);
        generator.setExcludeNodalPoint(nRootPoints > 0);
        generator.setMakeReeHover(nBody > 3);
        int startConnection = 1, stopConnection = 2;
        if (args.length > 2) {
            startConnection = Integer.parseInt(args[2]);
            stopConnection = Integer.parseInt(args[3]);
        }
        if (stopConnection > 0) {
            generator.setStartConnection(startConnection);
            generator.setStopConnection(stopConnection);
        }
        cluster.reset();
        generator.reset();
        long startTime = System.currentTimeMillis();
        int maxConnections = cluster.getNumConnections();
        int[] hist = new int[maxConnections];
        FileWriter fileWriter = null;
        String fileName = "";
        fileName = "B"+nBody;
        if (stopConnection > 0) {
            fileName += "_"+startConnection+"_"+stopConnection;
        }
        try {
            fileWriter = new FileWriter(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        if (startConnection == 0) {
            if (nBody > 3) {
                generator.calcReeHoover();
            }
            try {
                fileWriter.write(cluster.mReeHooverFactor+"\t"+cluster.toString()+"\n");
                fileWriter.flush();
            }
            catch(IOException e) {
                throw new RuntimeException("Cannot write to "+fileName+", caught IOException: "+e);
            }
            hist[maxConnections-1] = 1;
        }
        
        while (generator.advance()) {
            try {
                fileWriter.write(cluster.mReeHooverFactor+"\t"+cluster.toString()+"\n");
            }
            catch(IOException e) {
                throw new RuntimeException("Cannot write to "+fileName+", caught IOException: "+e);
            }
            hist[cluster.getNumConnections()-1]++;
        }
        try {
            fileWriter.close();
        }
        catch(IOException e) {
            throw new RuntimeException("Cannot close "+fileName+", caught IOException: "+e);
        }
        for (int i = 0; i < maxConnections; i++) {
            if (hist[i] > 0) {
                System.out.println(hist[i] + " clusters of size "+(i+1));
            }
        }
        System.out.println("total time " + (System.currentTimeMillis() - startTime) / 1000 + " seconds");
    }
    
    
    //work arrays used by methods that detect topology of cluster
    private final boolean[] mMarked;
    private final int[] mStack;
    
    //fields that flag types of clusters to be generated
    private boolean mOnlyDoublyConnected, mOnlyConnected;
    private boolean mExcludeNodalPoint, mExcludeArticulationPoint, mExcludeArticulationPair;
    private boolean mAllPermutations = true;
    private boolean excludeRootPermutations;
    private boolean mReeHoover;
    
    //used for Ree-Hoover cluster generation
    private ClusterDiagram mClusterCopy;
    
    
    private final int[] mConnectStack;
    private int mStartConnection, mStopConnection;
    
    //list of all potential connections in a cluster; effectively assigns index to each pair 
    private final int[][] mConnectList;
    private int mConnectStackLen;
    private int mTotalNumConnections;
    
    //cursor that holds the state between iterations
    private int mCurrentConnection;
    
    public final ClusterDiagram mCluster;
    private ClusterGenerator mReeHooverGenerator;


}
