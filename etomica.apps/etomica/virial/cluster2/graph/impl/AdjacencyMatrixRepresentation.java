package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.EdgesRepresentation;

/**
 * This class maps every edge of an N by N matrix onto a Bitmap using N^2 bits.
 * The edge (n1,n2) exists in the matrix iff the bit (n1*N + n2) is set. Since
 * this class assumes storage for both edges (n1,n2) and (n2,n1), it can be used
 * to represent edges for both directed and undirected graphs.
 * 
 * @author Demian Lessa
 */
public class AdjacencyMatrixRepresentation implements EdgesRepresentation {

  private byte nodeCount;

  public AdjacencyMatrixRepresentation(byte numNodes) {

    nodeCount = numNodes;
  }

  @Override
  public int getEdgeID(int fromNodeID, int toNodeID) {

    return fromNodeID * getNodeCount() + toNodeID;
  }

  @Override
  public int getFromNodeID(int edgeID) {

    return edgeID / getNodeCount();
  }

  @Override
  public int getToNodeID(int edgeID) {

    return edgeID % getNodeCount();
  }

  @Override
  public String toString(int edgeID) {

    return "(" + getFromNodeID(edgeID) + "," + getToNodeID(edgeID) + ")";
  }

  @Override
  public byte getNodeCount() {

    return nodeCount;
  }
}
