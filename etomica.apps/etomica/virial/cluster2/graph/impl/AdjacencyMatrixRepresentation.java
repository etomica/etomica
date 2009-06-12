package etomica.virial.cluster2.graph.impl;

/**
 * This class maps every edge of an N by N matrix onto a Bitmap using N^2 bits.
 * The edge (n1,n2) exists in the matrix iff the bit (n1*N + n2) is set. Since
 * this class assumes storage for both edges (n1,n2) and (n2,n1), it can be used
 * to represent edges for both directed and undirected graphs.
 * 
 * @author Demian Lessa
 */
public class AdjacencyMatrixRepresentation extends AbstractBitmapRepresentation {

  public AdjacencyMatrixRepresentation(byte numNodes) {

    super(numNodes);
  }

  public int getEdgeCount() {

    return getEdges().bitCount() / 2;
  }

  public int getEdgeID(int fromNodeID, int toNodeID) {

    if (fromNodeID > toNodeID) {
      return getEdgeID(toNodeID, fromNodeID);
    }
    return fromNodeID * getNodeCount() + toNodeID;
  }

  public int getFromNodeID(int edgeID) {

    return edgeID / getNodeCount();
  }

  public int getToNodeID(int edgeID) {

    return edgeID % getNodeCount();
  }
}