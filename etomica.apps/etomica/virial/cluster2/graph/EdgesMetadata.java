package etomica.virial.cluster2.graph;

public interface EdgesMetadata {

  public EdgesMetadata copy();

  public EdgesMetadata ncopy();

  public GraphCoefficient getCoefficient();

  public void setCoefficient(GraphCoefficient value);
}