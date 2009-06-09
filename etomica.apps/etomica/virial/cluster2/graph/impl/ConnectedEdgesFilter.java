package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.AbstractEdgesFilter;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.algorithms.GraphAlgorithmsFactory;
import etomica.virial.cluster2.graph.algorithms.GraphProperty;

public class ConnectedEdgesFilter extends AbstractEdgesFilter {

  private static final String FLAG_CONNECTED = "Connected";
  private GraphProperty connectedAlgorithm = GraphAlgorithmsFactory.getConnectedAlgo();
  private Nodes nodes;
  
  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************

  public ConnectedEdgesFilter(Nodes nodes) {

    super();
    this.nodes = nodes;
  }

  public ConnectedEdgesFilter(Nodes nodes, AbstractEdgesFilter filter) {

    super(filter);
    this.nodes = nodes;
  }

  // ***********************
  // * PROTECTED METHODS
  // ***********************

  @Override
  protected boolean doAccept(Edges edges) {

    return connectedAlgorithm.check(nodes, edges);
  }

  @Override
  protected boolean doPreAccept(Edges edges) {

    return true;
  }

  @Override
  protected String tag() {
    
    return ConnectedEdgesFilter.FLAG_CONNECTED;
  }
}