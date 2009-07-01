package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.EdgesMetadata;

public class SimpleEdgesMetadata implements EdgesMetadata {

  private double coefficient;

  public SimpleEdgesMetadata(double value) {

    coefficient = value;
  }

  public EdgesMetadata copy() {
    
    return new SimpleEdgesMetadata(coefficient);
  }
  
  public EdgesMetadata ncopy() {
    
    return new SimpleEdgesMetadata(-coefficient);
  }
  
  public double getCoefficient() {

    return coefficient;
  }

  public void setCoefficient(double value) {

    coefficient = value;
  }

  @Override
  public String toString() {

    return String.valueOf(Double.valueOf(coefficient).longValue());
  }
}