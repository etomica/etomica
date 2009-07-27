package etomica.virial.cluster2.graph;

public interface NodeAttributes {

  public char getColor();

  public boolean isCompatible(NodeAttributes attr);

  public boolean isSameClass(NodeAttributes attr);

  public boolean isSameColor(NodeAttributes attr);
}