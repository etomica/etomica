package etomica.graph.engine.impl;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.engine.Parser.Constructor;

public class ConstructorImpl implements Constructor {

  private Set<String> filters = new HashSet<String>();

  public void addFilter(String filterName) {

    this.filters.add(filterName);
  }

  public boolean hasFilter(String filterName) {

    return filters.contains(filterName);
  }
}