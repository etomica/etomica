package etomica.virial.cluster2.graph;

public interface Nodes {

  // how many nodes are there?
  public byte count();
  
  // get the attributes of nodeID
  public NodeAttributes getAttributes(int nodeID);
}