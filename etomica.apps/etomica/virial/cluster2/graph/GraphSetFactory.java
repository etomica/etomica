package etomica.virial.cluster2.graph;

import etomica.virial.cluster2.graph.impl.SimpleNodes;


public class GraphSetFactory {
  
  public static GraphSet completeGraphSet(char[] fieldColors, char[] rootColors) {

    Nodes nodes = new SimpleNodes(fieldColors, rootColors);
    return GraphFactory.storedGraphSet(nodes, null, false);

  }

  public static GraphSet isomorphFreeGraphSet(char[] fieldColors, char[] rootColors) {

    Nodes nodes = new SimpleNodes(fieldColors, rootColors);
    return GraphFactory.storedGraphSet(nodes, null, true);
  }
}