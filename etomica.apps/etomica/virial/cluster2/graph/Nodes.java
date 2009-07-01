package etomica.virial.cluster2.graph;

import java.util.List;

public interface Nodes {

  // how many nodes are there?
  public byte count();

  // get the attributes of nodeID
  public NodeAttributes getAttributes(int nodeID);

  public List<Integer> getPartition(int partitionID);

  public int getPartitionCount();
}