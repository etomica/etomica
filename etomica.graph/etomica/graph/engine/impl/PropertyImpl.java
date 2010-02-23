package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Property;

public class PropertyImpl implements Property {

  private String property;

  public PropertyImpl(String property) {

    this.property = property;
  }

  public String getName() {

    return this.property;
  }
}
