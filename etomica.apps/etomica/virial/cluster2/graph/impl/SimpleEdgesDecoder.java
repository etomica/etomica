package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.EdgesDecoder;

/*
 * This decoder maps to and from a set of edges encoded as a bitmap. It assumes
 * the identity mapping from bitIndex to edgeID, resulting in the following
 * mapping between edges and bits in the bitmap:
 * 
 * bitIndex  edgeID  (fromNodeID,toNodeID)
 * --------  ------  ---------------------
 *    0         0      (0,1)
 *    1         1      (0,2)
 *    2         2      (0,3)
 *    3         3      (1,2)
 *    4         4      (1,3)
 *    5         5      (2,3)
 * 
 * @author Demian Lessa
 * 
 */
public class SimpleEdgesDecoder implements EdgesDecoder {

  private byte nodeCount;

  public SimpleEdgesDecoder(byte numNodes) {

    nodeCount = numNodes;
  }

  /**
   * Returns the actual number of bits used to represent this set of edges.
   * 
   */

  @Override
  public int getBitCapacity() {

    return nodeCount * (nodeCount - 1) / 2;
  }

  @Override
  public int getBitIndex(int edgeID) {

    return edgeID;
  }

  @Override
  public int getBitIndex(int fromNodeID, int toNodeID) {

    return getBitIndex(getEdgeID(fromNodeID, toNodeID));
  }

  @Override
  public int getEdgeID(int bitIndex) {

    return bitIndex;
  }

  /*
   * bitIndex  edgeID  (fromNodeID,toNodeID)
   * --------  ------  ---------------------
   *    0         0      (0,1)
   *    1         1      (0,2)
   *    2         2      (0,3)
   *    3         3      (1,2)
   *    4         4      (1,3)
   *    5         5      (2,3)
   * (non-Javadoc)
   * @see etomica.virial.cluster2.graph.EdgesDecoder#getEdgeID(int, int)
   */
  @Override
  public int getEdgeID(int fromNodeID, int toNodeID) {

    if (fromNodeID > toNodeID) {
      return getEdgeID(toNodeID, fromNodeID);
    }
    return (toNodeID - fromNodeID - 1) + sumMaxEdges(0, fromNodeID - 1);
  }

  @Override
  public int getFromNodeID(int bitIndex) {

    int fromNodeID = 0;
    int offset = maxEdges(fromNodeID) - 1;
    int edgeID = getEdgeID(bitIndex);
    while (edgeID > offset) {
      fromNodeID++;
      offset += maxEdges(fromNodeID);
    }
    return fromNodeID;
  }

  @Override
  public int getToNodeID(int bitIndex) {

    int edgeID = getEdgeID(bitIndex);
    int fromNodeID = getFromNodeID(bitIndex);
    int fromEdgeID = getEdgeID(fromNodeID, fromNodeID + 1);
    return (fromNodeID + 1) + (edgeID - fromEdgeID);
  }

  @Override
  public String toString(int bitIndex) {

    return "(" + getFromNodeID(bitIndex) + "," + getToNodeID(bitIndex) + ")";
  }

  /*
   * Returns the largest number of edges encoded with the fromNodeID as the
   * first node in the edge pair (fromNodeID, toNodeID). For N=4 we have the
   * following:
   * 
   * fromNodeID maxEdges maxEdgeSet
   * ---------- -------- --------------------- 
   *         0        3 <(0,1), (0,2), (0,3)> 
   *         1        2 <(1,2), (1,3)> 
   *         2        1 <(2,3)> 
   *         3        0
   * 
   */
  protected int maxEdges(int fromNodeID) {

    return (nodeCount - fromNodeID - 1);
  }

  protected byte nodeCount() {

    return nodeCount;
  }

  protected int sumMaxEdges(int firstNodeID, int lastNodeID) {

    int result = 0;
    for (int nodeID = firstNodeID; nodeID <= lastNodeID; nodeID++) {
      result += maxEdges(nodeID);
    }
    return result;
  }
}