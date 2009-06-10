package etomica.virial.cluster2.graph;

/*
 * This interface provides a simple mechanism to map pairs of nodes onto a
 * linear scale of edge identifiers. Hence, it provides a means to access
 * edge information from some particular storage that is indexed on a single
 * integer, such as a vector or a bitmap.
 * 
 * @author Demian Lessa
 * 
 */
public interface EdgesRepresentation {

  public int getEdgeID(int fromNodeID, int toNodeID);

  public int getFromNodeID(int edgeID);

  public byte getNodeCount();

  public int getToNodeID(int edgeID);

  public String toString(int edgeID);
}