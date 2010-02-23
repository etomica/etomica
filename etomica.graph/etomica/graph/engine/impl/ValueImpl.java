package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Value;

public class ValueImpl implements Value {

  private String value;

  public ValueImpl(String value) {

    this.value = value;
  }

  public String getValue() {

    return this.value;
  }
}