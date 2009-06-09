package etomica.virial.cluster2.graph.algorithms;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.NodesVisitor;

public class GraphAlgorithms {

  public static boolean hasArticulationPair(final Nodes nodes, final Edges edges) {

    return GraphAlgorithmsFactory.getArticulationPairAlgo().check(nodes, edges);
  }

  public static boolean hasArticulationPoint(final Nodes nodes,
      final Edges edges) {

    return GraphAlgorithmsFactory.getArticulationPointAlgo()
        .check(nodes, edges);
  }

  public static boolean hasNodalPoint(final Nodes nodes, final Edges edges) {

    return GraphAlgorithmsFactory.getNodalPointAlgo().check(nodes, edges);
  }

  public static boolean isBiconnected(final Nodes nodes, final Edges edges) {

    return GraphAlgorithmsFactory.getBiconnectedAlgo().check(nodes, edges);
  }

  public static boolean isConnected(final Nodes nodes, final Edges edges) {

    return GraphAlgorithmsFactory.getConnectedAlgo().check(nodes, edges);
  }

  public static boolean checkIsomorphic(final Nodes nodes, final Edges edges1,
      final Edges edges2) {

    return GraphAlgorithmsFactory.getIsomorphismAlgo().check(nodes, edges1,
        edges2);
  }

  public static boolean traverseBF(int nodeID, final Nodes nodes,
      final Edges edges, final NodesVisitor visitor) {

    return GraphAlgorithmsFactory.getBFTAlgo().traverseComponent(nodeID, nodes,
        edges, visitor);
  }

  public static void traverseBF(final Nodes nodes, final Edges edges,
      final NodesVisitor visitor) {

    GraphAlgorithmsFactory.getBFTAlgo().traverseAll(nodes, edges, visitor);
  }

  public static boolean traverseDF(int nodeID, final Nodes nodes,
      final Edges edges, final NodesVisitor visitor) {

    return GraphAlgorithmsFactory.getDFTAlgo().traverseComponent(nodeID, nodes,
        edges, visitor);
  }

  public static void traverseDF(final Nodes nodes, final Edges edges,
      final NodesVisitor visitor) {

    GraphAlgorithmsFactory.getDFTAlgo().traverseAll(nodes, edges, visitor);
  }
}