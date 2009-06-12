package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.EdgesMetadata;

public class SimpleEdgesMetadata implements EdgesMetadata {

  private double coefficient;

  public SimpleEdgesMetadata(double value) {

    coefficient = value;
  }

  public double getCoefficient() {

    return coefficient;
  }

  public void setCoefficient(double value) {

    coefficient = value;
  }
}