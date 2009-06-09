package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.EdgesMetadata;

public class SimpleEdgesMetadata implements EdgesMetadata {

  private double coefficient;

  public SimpleEdgesMetadata(double value) {

    coefficient = value;
  }

  @Override
  public double getCoefficient() {

    return coefficient;
  }

  @Override
  public void setCoefficient(double value) {

    coefficient = value;
  }
}