package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.EdgesMetadata;
import etomica.virial.cluster2.graph.GraphCoefficient;

public class SimpleEdgesMetadata implements EdgesMetadata {

  private GraphCoefficient coefficient;

  public SimpleEdgesMetadata(GraphCoefficient value) {

    coefficient = value;
  }

  public EdgesMetadata copy() {
    
    return new SimpleEdgesMetadata(coefficient);
  }
  
  public EdgesMetadata ncopy() {
    
    return new SimpleEdgesMetadata(coefficient.switchSign());
  }
  
  public GraphCoefficient getCoefficient() {

    return coefficient;
  }

  public void setCoefficient(GraphCoefficient value) {

    coefficient = value;
  }

  @Override
  public String toString() {

    return coefficient.toString();
  }
}