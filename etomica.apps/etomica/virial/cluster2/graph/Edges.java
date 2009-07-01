package etomica.virial.cluster2.graph;

public interface Edges {

  // returns the canonically labeled set of edges isomorphic to this set of edges
  public Edges canonical();

  // returns the set of edges corresponding to the complement of this set of edges
  public Edges complement();

  // returns an independent copy of this set of edges
  public Edges copy();

  // returns an independent copy of this set of edges with the inverse of the coefficient
  public Edges ncopy();

  // returns the number of edges
  public int count();

  // get the attributes of the edge (node1, node2)
  public EdgeAttributes getAttributes(int fromNodeID, int toNodeID);

  // returns the number of head nodes adjacent to nodeID
  public int getInDegree(int nodeID);

  // returns the i-th head node adjacent to nodeID
  public int getInNode(int nodeID, int index);

  public EdgesMetadata getMetadata();

  // returns the number of tail nodes adjacent to nodeID
  public int getOutDegree(int nodeID);

  // returns the i-th tail node adjacent to nodeID
  public int getOutNode(int nodeID, int index);

  public boolean hasEdge(int fromNodeID, int toNodeID);

  // returns the object with the internal representation of the edges
  public EdgesRepresentation getRepresentation();
}