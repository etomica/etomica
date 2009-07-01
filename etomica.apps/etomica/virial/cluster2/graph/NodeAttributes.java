package etomica.virial.cluster2.graph;

public interface NodeAttributes {

  public boolean isCompatible(NodeAttributes attr);

  public boolean isSameColor(NodeAttributes attr);
}