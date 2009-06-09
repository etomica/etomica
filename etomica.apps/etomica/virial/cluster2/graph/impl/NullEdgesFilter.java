package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.AbstractEdgesFilter;
import etomica.virial.cluster2.graph.Edges;

public class NullEdgesFilter extends AbstractEdgesFilter {

  private static final String FLAG_NULL = "NULL";

  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************

  public NullEdgesFilter() {

    super();
  }

  public NullEdgesFilter(AbstractEdgesFilter filter) {

    super(filter);
  }

  // ***********************
  // * PROTECTED METHODS
  // ***********************

  @Override
  protected boolean doAccept(Edges edges) {

    return false;
  }

  @Override
  protected boolean doPreAccept(Edges edges) {

    return false;
  }

  @Override
  protected String tag() {
    
    return NullEdgesFilter.FLAG_NULL;
  }
}