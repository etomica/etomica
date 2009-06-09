package etomica.virial.cluster2.graph;

/*
 * This interface decodes the internal representation of the edges of a graph
 * into a consistent form across implementations. This adds immediate support
 * for arbitrary edges generators. Each generator must implement a decoder to
 * translate its internal representation to the standard one.
 * 
 * The SimpleEdgesDecoder implements this interface and encodes edges in the
 * range 0..N(N-1)/2-1. Edges are represented in a packed, upper triangular 
 * adjacency matrix. For instance, for N=4, a possible mapping between edgeID 
 * and edge coordinates is given below:
 * 
 * edgeID  (fromNodeID,toNodeID)
 * ------  ---------------------
 *    0      (0,1)
 *    1      (0,2)
 *    2      (0,3)
 *    3      (1,2)
 *    4      (1,3)
 *    5      (2,3)
 * 
 * @author Demian Lessa
 * 
 */
public interface EdgesDecoder {

  public int getBitCapacity();

  public int getBitIndex(int edgeID);

  public int getBitIndex(int fromNodeID, int toNodeID);

  public int getEdgeID(int bitIndex);

  public int getEdgeID(int fromNodeID, int toNodeID);

  public int getFromNodeID(int bitIndex);

  public int getToNodeID(int bitIndex);
  
  public String toString(int bitIndex);
}