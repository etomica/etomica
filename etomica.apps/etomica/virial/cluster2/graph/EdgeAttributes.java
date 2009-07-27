package etomica.virial.cluster2.graph;

public interface EdgeAttributes {

  public char getColor();

  public boolean isCompatible(EdgeAttributes attr);

  public boolean isSameColor(EdgeAttributes attr);
}
