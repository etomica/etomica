package etomica.modules.clustergenerator;



/**
 * @author andrew
 */
public class GenCluster {
    public final int mNumBody;
    private final boolean[] mMarked;
    public final int[][] mConnections;
    private final int[] mStack;
    public final int[] mScore;
    private GenCluster mCopy;
    public int mNumIdenticalPermutations;
    private boolean mOnlyDoublyConnected, mOnlyConnected;
    private boolean mExcludeNodalPoint, mExcludeArticulationPoint, mExcludeArticulationPair;
    private boolean mAllPermutations = true;
    private int[] mConnectStack;
    private int[][] connectList;
    private int mConnectStackLen;
    private int mCurrentConnection;
    private int mNumConnections;
    private int mNumRootPoints;
    private int[] mRootPoints;
    private final boolean[] mIsRootPoint;

    public GenCluster(int aNBody) {
        mNumBody = aNBody;
        mMarked = new boolean[mNumBody];
        mConnections = new int[mNumBody][mNumBody];
        makeFullStar();
        mStack = new int[mNumBody];
        mIsRootPoint = new boolean[mNumBody];
        setNumRootPoints(0);
        mScore = new int[mNumBody];
    }

    public GenCluster(GenCluster cluster) {
        mNumBody = cluster.mNumBody;
        mMarked = new boolean[mNumBody];
        mConnections = new int[mNumBody][];
        mIsRootPoint = new boolean[mNumBody];
        mStack = new int[mNumBody];
        mRootPoints = new int[mNumRootPoints];
        for (int i = 0; i < mNumBody; i++) {
            mConnections[i] = new int[mNumBody];
            System.arraycopy(cluster.mConnections[i], 0, mConnections[i], 0,
                    mNumBody);
        }
        System.arraycopy(cluster.mIsRootPoint, 0, mIsRootPoint, 0, mNumBody);
        System.arraycopy(cluster.mRootPoints, 0, mRootPoints, 0, mNumRootPoints);
        mScore = new int[mNumBody];
        System.arraycopy(cluster.mScore, 0, mScore, 0, mNumBody);
        mNumIdenticalPermutations = cluster.mNumIdenticalPermutations;
    }

    public void setOnlyConnected(boolean flag) {
        mOnlyConnected = flag;
    }

    public void setOnlyDoublyConnected(boolean flag) {
        mOnlyDoublyConnected = flag;
    }

    public void setExcludeNodalPoint(boolean flag) {
        mExcludeNodalPoint = flag;
    }

    public void setExcludeArticulationPoint(boolean flag) {
        mExcludeArticulationPoint = flag;
    }

    public void setExcludeArticulationPair(boolean flag) {
        mExcludeArticulationPair = flag;
    }

    public void setAllPermutations(boolean flag) {
        mAllPermutations = flag;
        if (!mAllPermutations) {
            mCopy = new GenCluster(mNumBody);
        }
    }

    public void setNumRootPoints(int n) {
        mNumRootPoints = n;
        mRootPoints = new int[mNumRootPoints];
        for (int i=0; i<mNumRootPoints; i++) {
            mRootPoints[i] = i;
            mIsRootPoint[i] = true;
        }
    }
    public int getNumRootPoints() {
        return mNumRootPoints;
    }
    
    /**
     * returns the number of active connections in the cluster.  
     * assumes connections were deleted via advance()
     * @return
     */
    public int getNumConnections() {
        return mNumConnections - mConnectStackLen;
    }
    
    /**
     * turns this cluster into a full star -- each point in the cluster is
     * connected to every other point in the cluster.
     */
    public void makeFullStar() {
        int i, j, k;
        for (i = 0; i < mNumBody; i++) {
            k = 0;
            for (j = 0; j < mNumBody; j++) {
                if (i != j) {
                    mConnections[i][k] = j;
                    k++;
                }
            }
            mConnections[i][mNumBody - 1] = -1;
        }
        mNumConnections = mNumBody*(mNumBody-1)/2;
    }

    /**
     * Copies the only connections to another cluster.
     * 
     * @param dst cluster to receive the connections
     */
    public void copy(GenCluster dst) {
        for (int i = 0; i < mNumBody; i++) {
            System.arraycopy(mConnections[i], 0, dst.mConnections[i], 0, mNumBody - 1);
        }
        dst.setNumRootPoints(mNumRootPoints);
        dst.mNumRootPoints = mNumRootPoints;
        dst.mNumIdenticalPermutations = mNumIdenticalPermutations;
    }

    public void swap(int point1, int point2) {
        // first have things pointing at point1 to point at point2
        // and vica-versa
        for (int i = 0; i < mNumBody; i++) {
            boolean foundPoint1, foundPoint2;
            foundPoint1 = false;
            foundPoint2 = false;
            for (int j = 0; mConnections[i][j] != -1; j++) {
                if (!foundPoint1 && mConnections[i][j] == point1) {
                    mConnections[i][j] = point2;
                    foundPoint1 = true;
                    if (foundPoint2) {
                        break;
                    }
                }
                else if (!foundPoint2 && mConnections[i][j] == point2) {
                    mConnections[i][j] = point1;
                    foundPoint2 = true;
                    if (foundPoint1) {
                        break;
                    }
                }
            }
        }
        // actually swap the list of connections for point1 and point2
        int[] tmpConnections = new int[mNumBody - 1];
        System.arraycopy(mConnections[point1], 0, tmpConnections, 0, mNumBody - 1);
        System.arraycopy(mConnections[point2], 0, mConnections[point1], 0, mNumBody - 1);
        System.arraycopy(tmpConnections, 0, mConnections[point2], 0, mNumBody - 1);
    }

    /**
     * Determines whether any cluster permuted from the current one by
     * exchanging points from startPoint to stopPoint (inclusive) has a higher
     * score than cluster2. If such a cluster is located, the method returns
     * true immeadiately. The method always returns the connections to their
     * state when the method was called.  If earlyReturn is true, the method
     * returns false if any of the permuted clusters is identical to cluster2.
     * 
     * @param startPoint starting point for range over which to perform permutations
     * @param stopPoint stopping point for ranger over which to perform permutations
     * @param cluster2 cluster to compare each permuted cluster with
     * @return true if any permuted cluster has a higher score than cluster2
     */
    private boolean swapAllAndCompare(int startPoint, int stopPoint, GenCluster cluster2) {
        // recursively swap everything between startPoint+1 and stopPoint
        // this is basically the 0th iteration below
        if (startPoint+1 < stopPoint && swapAllAndCompare(startPoint + 1, stopPoint, cluster2)) {
            return true;
        }
        // recursively swap startPoint with startPoint+1 to stopPoint
        for (int i = startPoint + 1; i < stopPoint + 1; i++) {
            if (mIsRootPoint[startPoint] == mIsRootPoint[i]) {
                swap(startPoint, i);
                boolean isGreaterThan = scoreGreaterThan(cluster2);
                isGreaterThan = isGreaterThan || swapAllAndCompare(startPoint + 1, stopPoint, cluster2);
                swap(startPoint, i);
                if (isGreaterThan) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Prints all permutations of the current cluster permuting points between
     * the given starting and stopping points.
     * 
     * @param startPoint starting point for range over which to perform permutations
     * @param stopPoint stopping point for ranger over which to perform permutations
     */
    private void printAllPermutations(int startPoint, int stopPoint) {
        if (stopPoint > startPoint + 1) {
            printAllPermutations(startPoint + 1, stopPoint);
            for (int i = startPoint + 1; i < stopPoint + 1; i++) {
                swap(startPoint, i);
                printAllPermutations(startPoint + 1, stopPoint);
                swap(startPoint, i);
            }
        }
        else {
            System.out.println(toString());
            swap(startPoint, stopPoint);
            System.out.println(toString());
            swap(startPoint, stopPoint);
        }
    }

    /**
     * determines if the cluster is connected after removal of a node
     * 
     * @param exclude node to be removed
     * @return true if the cluster is still connected after node removal
     */
    public boolean isConnected(int exclude) {
        for (int i = 0; i < mNumBody; i++) {
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
     * determines if the cluster is connected
     */
    public boolean isConnected() {
        for (int i = 0; i < mNumBody; i++) {
            mMarked[i] = false;
        }
        mMarked[0] = true;
        mStack[0] = 0;
        mStack[1] = -1;
        return isConnectedReally(1);
    }
    
    private boolean isConnectedReally(int numMarked) {
        int stackLength = 1;
        while (stackLength > 0) {
            // pop the mStack
            int thisNode = mStack[stackLength - 1];
            stackLength--;
            for (int i = 0;; i++) {
                int connectedNode = mConnections[thisNode][i];
                if (connectedNode == -1) {
                    // end of connections for thisNode
                    break;
                }
                if (!mMarked[connectedNode]) {
                    mMarked[connectedNode] = true;
                    numMarked++;
                    if (numMarked == mNumBody) {
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
     * convenience method to determine whether the cluster is elementary or not
     * 
     * @return true if the cluster is elementary
     */
    public boolean isElementary() {
        for (int i = mNumBody - 1; i >= 0; i--) {
            if (isNodalPoint(i) || isArticulationPoint(i)) return false;
            for (int j = i - 1; j >= 0; j--) {
                if (isArticulationPair(i, j)) return false;
            }
        }
        return true;
    }

    /**
     * determines if the given node is a nodal point
     * 
     * @return true if the point is a nodal point
     */
    public boolean isNodalPoint(int exclude) {
        if (mIsRootPoint[exclude] || mNumRootPoints < 2) {
            return false;
        }
        int stackLength = 1;
        int i;
        int thisNode, connectedNode;
        int rootPointsFound = 1;

        for (i = 0; i < mNumBody; i++) {
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
                connectedNode = mConnections[thisNode][i];
                if (connectedNode == -1) {
                    // end of connections for thisNode
                    break;
                }
                if (!mMarked[connectedNode]) {
                    if (mIsRootPoint[connectedNode]) {
                        if (++rootPointsFound == mNumRootPoints) {
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
     * determines if the given point is an articulation point
     * 
     * @param exclude point to be checked as an articulation point
     * @return true if the passed argument is an articulation point
     */
    public boolean isArticulationPoint(int exclude) {
        int i;
        int numMarksNeeded;
        // # of root points that haven't been excluded
        int nRootPoints = mNumRootPoints;;

        for (i = 0; i < mNumBody; i++) {
            mMarked[i] = false;
        }
        mMarked[exclude] = true;
        numMarksNeeded = mNumBody - mNumRootPoints;
        if (!mIsRootPoint[exclude]) {
            numMarksNeeded--;
        }
        else {
            // excluded a root point, so don't make isArticulated try to find it
            nRootPoints--;
        }
        mStack[0] = mNumRootPoints;
        while (mMarked[mStack[0]]) {
            mStack[0]++;
        }
        mMarked[mStack[0]] = true;
        mStack[1] = -1; // -1 terminates the mStack

        return isArticulated(numMarksNeeded,nRootPoints);
    }

    /**
     * determines if the given pair of points is an articulation pair
     * 
     * @param exclude1 first point of pair
     * @param exclude2 second point of pair
     * @return true if the pair is an articulation point
     */
    public boolean isArticulationPair(int exclude1, int exclude2) {
        int i;
        int numMarksNeeded;
        // # of root points that haven't been excluded
        int nRootPoints = mNumRootPoints;

        for (i = 0; i < mNumBody; i++) {
            mMarked[i] = false;
        }
        // numMarked only includes marked filled points
        numMarksNeeded = mNumBody - mNumRootPoints;
        mMarked[exclude1] = true;
        if (!mIsRootPoint[exclude1]) {
            numMarksNeeded--;
        }
        else {
            // excluded a root point, so don't make isArticulated try to find it
            nRootPoints--;
        }
        mMarked[exclude2] = true;
        if (!mIsRootPoint[exclude2]) {
            numMarksNeeded--;
        }
        else {
            nRootPoints--;
        }
        mStack[0] = mNumRootPoints;
        while (mMarked[mStack[0]]) {
            mStack[0]++;
            if (mStack[0] == mNumBody) {
                // excluded all filled points
                return false;
            }
        }
        mMarked[mStack[0]] = true;
        mStack[1] = -1; // -1 terminates the mStack
        
        return isArticulated(numMarksNeeded, nRootPoints);
    }

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
                    int connectedNode = mConnections[thisNode][i];
                    if (connectedNode == -1) {
                        // end of connections for thisNode
                        break;
                    }
                    if (!mMarked[connectedNode]) {
                        if (!mIsRootPoint[connectedNode]) {
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
            int thisNode = mNumRootPoints+1;
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
     * deletes a connection from one node of the cluster to another. the method
     * should be called again to delete the connection from the second node to
     * the first
     * 
     * @param node1 node owning the connection to be deleted
     * @param node2 node connected to node1
     * @return true if the deletion was successful.
     */
    public boolean deleteConnection(int node1, int node2) {
        boolean success = false;
        int[] connections1 = mConnections[node1];
        for (int i = 0; ; i++) {
            int connection = connections1[i];
            if (connection == -1) {
                // ran out of connections. we're done
                break;
            }
            // did we find it?
            success = success || connection == node2;
            if (success) {
                // shift
                mConnections[node1][i] = mConnections[node1][i + 1];
            }
        }
        return success;
    }

    /**
     * adds a connection between the given nodes of the cluster.
     * The method should be called again to create the connection
     * from the second node to the first 
     * @param node1 node owning the connection to be added
     * @param node2 node to be connected to node1
     */
    public void addConnection(int node1, int node2) {
        int[] connections1 = mConnections[node1];
        for (int i = 0;; i++) {
            if (connections1[i] == node2) {
                throw new RuntimeException("attempted to add already-existing connection");
            }
            if (connections1[i] == -1) {
                connections1[i] = node2;
                return;
            }
        }
    }

    public String toString() {
        String out = "";
        int i, j;
        for (i = 0; i < mNumBody; i++) {
            for (j = 0; j < mNumBody - 1 && mConnections[i][j] != -1; j++) {
                if (i < mConnections[i][j]) {
                    out += i + " " + mConnections[i][j] + "\t";
                }
            }
        }
        return out;
    }

    /**
     * Prints all permutations of the current cluster.
     */
    public void printAllPermutations() {
        if (mCopy == null) {
            mCopy = new GenCluster(mNumBody);
        }
        copy(mCopy);
        mCopy.printAllPermutations(1, mNumBody - 2);
        /*
         * copy.swap(0,mNumBody-1); copy.printAllPermutations(1,mNumBody-2);
         * copy.swap(0,mNumBody-1);
         */
    }

    /**
     * Determines if the cluster has a score less than that of another cluster.
     * The score of the other cluster is computed only enough to determine if it
     * is less than the current cluster.
     * 
     * @return true if the cluster's score is less than the given one.
     */
    private boolean scoreGreaterThan(GenCluster cluster2) {
        for (int i = 1; i < mNumBody / 2 + 1; i++) {
            int myScore = 0;
            for (int thisNode = 0; thisNode < mNumBody; thisNode++) {
                int thisScore = 1 << (mNumBody - 1 - thisNode);
                int target = thisNode + i;
                if (target > mNumBody - 1) target -= mNumBody;
                for (int j = 0; mConnections[thisNode][j] != -1; j++) {
                    if (mConnections[thisNode][j] == target) {
                        myScore += thisScore;
                        if (myScore > cluster2.mScore[i]) {
                            // the given cluster has a lower score at level i
                            return true;
                        }
                        break;
                    }
                }
            }
            if (myScore < cluster2.mScore[i]) {
                return false;
            }
        }
        cluster2.mNumIdenticalPermutations++;
        return false;
    }

    /**
     * compute a score for the cluster. The score is returned as an array of
     * integers. The first element is the score for connections between i and
     * i+1 and the second element is the score for connections between i and i+2
     * (the zeroth element is always 0). The score puts more importance on early
     * points (having a connection for point 0 will score one higher than having
     * the same connection for every other point).
     * 
     * @return cluster score as an array of scores for each type of connection
     */
    private void calcScore() {
        for (int thisNode = 1; thisNode < 1 + mNumBody/2; thisNode++) {
            mScore[thisNode] = 0;
        }
        for (int thisNode = 0; thisNode < mNumBody; thisNode++) {
            int thisScore = 1 << (mNumBody - 1 - thisNode);
            for (int i = 1; i < mNumBody / 2 + 1; i++) {
                int target = thisNode + i;
                if (target > mNumBody - 1) target -= mNumBody;
                for (int j = 0; mConnections[thisNode][j] != -1; j++) {
                    if (mConnections[thisNode][j] == target) {
                        mScore[i] += thisScore;
                        break;
                    }
                }
            }
        }
    }

    /**
     * Determines if the cluster has a maximum score. If the score
     * can be increased by swapping filled points or by swapping reference
     * nodes, then one of its permutations must have a higher score and the
     * method returns false
     * 
     * @return true if the cluster has the maximum score.
     */
    /*
     * this method assumes that an elementary cluster with a maximum score will
     * have score(node_0) >= score(node_mNumBody-1) and score(node_i) >=
     * score(node_i+1) for 0 < i < mNumBody-1
     */
    public boolean isMaximumScore() {
        calcScore();
        copy(mCopy);
        mNumIdenticalPermutations = 1;
/*        int[] pointScore = new int[mNumBody];
        for (int i=0; i<mNumBody; i++) {
            mCopy.mIdenticalTo[i] = i;
        }
        for (int i=0; i<mNumBody; i++) {
            int iScore = 0;
            int[] iConnections = mConnections[i];
            for (int j=0; iConnections[j] != -1; j++) {
                iScore += 1 << iConnections[j];
            }
            pointScore[i] = iScore;
            for (int j=0; j<i; j++) {
                if ((iScore == pointScore[j] || iScore-(1<<j) == pointScore[j]-(1<<i))
                        && mIsRootPoint[i] == mIsRootPoint[j]) {
                    mCopy.mIdenticalTo[i] = j;
                    int numIdenticalPoints = 0;
                    // find how many points identical to i
                    for (int k=0; k<mNumBody; k++) {
                        if (mCopy.mIdenticalTo[k] == mCopy.mIdenticalTo[i]) {
                            numIdenticalPoints++;
                        }
                    }
                    mNumIdenticalPermutations *= numIdenticalPoints;
                    break;
                }
            }
        }*/
        return !mCopy.swapAllAndCompare(0, mNumBody - 1, this);
    }

    public void reset() {
        makeFullStar();

        for (int i=mNumRootPoints-1; i>0; i--) {
            mNumConnections -= i;
        }
        for (int i=0; i<mNumRootPoints; i++) {
            for (int j=0; j<i; j++) {
                deleteConnection(i,j);
                deleteConnection(j,i);
            }
        }

        connectList = new int[mNumConnections][2];
        for (int k = 0, i = 0; i < mNumBody-1; i++) {
            for (int j = i + 1; j < mNumBody; j++) {
                if (!mIsRootPoint[i] || !mIsRootPoint[j]) {
                    connectList[k][0] = i;
                    connectList[k][1] = j;
                    k++;
                }
            }
        }

        mConnectStack = new int[mNumConnections];
        mCurrentConnection = 0;
        mConnectStackLen = 0;
        if (!mAllPermutations && !isMaximumScore()) {
            throw new RuntimeException("initial cluster not maxed");
        }
    }

    public boolean advance() {
        for (; mCurrentConnection < mNumConnections || mConnectStackLen > 0; mCurrentConnection++) {
            boolean success = true;
            if (mCurrentConnection == mNumConnections) {
                // we made it through the whole list. pop a connection off the
                // mStack and go through the list again from there.
                mCurrentConnection = mConnectStack[mConnectStackLen - 1];
                int node1 = connectList[mCurrentConnection][0];
                int node2 = connectList[mCurrentConnection][1];
                // System.out.println("adding back connection between "+node1+" and "+node2);
                mConnectStackLen--;
                addConnection(node1, node2);
                addConnection(node2, node1);
                continue;
            }
            for (int l = 0; l < mConnectStackLen; l++) {
                if (mConnectStack[l] == mCurrentConnection) {
                    // the connection has already been deleted
                    continue;
                }
            }
            int node1 = connectList[mCurrentConnection][0];
            int node2 = connectList[mCurrentConnection][1];
            // System.out.println("deleting connection betweeen "+node1+" and "+node2);
            if (!deleteConnection(node1, node2)) {
                throw new RuntimeException(
                        "failed to delete a connection that should have existed");
            }
            if (!deleteConnection(node2, node1)) {
                throw new RuntimeException(
                        "failed to delete a connection that should have existed");
            }
            if (mExcludeNodalPoint) {
                for (int i = mNumBody - 1; i >= 0; i--) {
                    if (isNodalPoint(i)) {
                        success = false;
                        break;
                    }
                }
            }
            if (success && mExcludeArticulationPoint) {
                for (int i = mNumBody - 1; i >= 0; i--) {
                    if (isArticulationPoint(i)) {
                        success = false;
                        break;
                    }
                }
            }
            if (success && mExcludeArticulationPair) {
                for (int i = mNumBody - 1; i >= 0; i--) {
                    for (int j = i - 1; j >= 0; j--) {
                        if (isArticulationPair(i, j)) {
                            success = false;
                            break;
                        }
                    }
                }
            }
            if (success && mOnlyDoublyConnected) {
                for (int i = mNumBody - 1; i >= 0; i--) {
                    if (!isConnected(i)) {
                        success = false;
                        break;
                    }
                }
            }
            if (success && mOnlyConnected) {
                success = isConnected();
            }
            if (success) {
                mConnectStack[mConnectStackLen] = mCurrentConnection;
                mConnectStackLen++;
                if (mAllPermutations || isMaximumScore()) {
                    mCurrentConnection++;
                    return true;
                }
            }
            else {
                addConnection(node1, node2);
                addConnection(node2, node1);
            }
        }
        return false;
    }


    public static void main(String[] args) {
        int nBody = 8;
        GenCluster cluster = new GenCluster(nBody);
        cluster.setAllPermutations(false);
        cluster.setOnlyDoublyConnected(false);
        cluster.setExcludeArticulationPoint(true);
        cluster.setExcludeArticulationPair(true);
        cluster.setExcludeNodalPoint(true);
        cluster.setNumRootPoints(2);
        cluster.reset();
        long startTime = System.currentTimeMillis();
        int maxConnections = cluster.getNumConnections();
        int[] hist = new int[maxConnections];
        hist[maxConnections-1] = 1;
        while (cluster.advance()) {
//            System.out.println(cluster.mNumIdenticalPermuations + "x\t"+cluster.toString());
            hist[cluster.getNumConnections()-1]++;
        }
        for (int i = 0; i < maxConnections; i++) {
            if (hist[i] > 0) {
                System.out.println(hist[i] + " clusters of size "+i);
            }
        }
        System.out.println("total time " + (System.currentTimeMillis() - startTime) / 1000 + " seconds");
    }

}