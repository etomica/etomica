package etomica.virial.cluster2.graph;

public interface EdgesMetadata {

  public EdgesMetadata copy();

  public EdgesMetadata ncopy();

  public double getCoefficient();

  public void setCoefficient(double value);
}