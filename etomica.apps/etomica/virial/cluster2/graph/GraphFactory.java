package etomica.virial.cluster2.graph;

import java.util.ArrayList;
import java.util.List;

import etomica.virial.cluster2.graph.impl.NaiveEdgesGenerator;
import etomica.virial.cluster2.graph.impl.NautyEdgesGenerator;
import etomica.virial.cluster2.graph.impl.SimpleEdges;
import etomica.virial.cluster2.graph.impl.SimpleEdgesMetadata;
import etomica.virial.cluster2.graph.impl.SimpleGraphSet;
import etomica.virial.cluster2.nauty.NautyInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;

public class GraphFactory {

  // default representation is adjacency matrix
  public static boolean USE_UPPER_TRIANGLE = true;
  // a class of edge attributes that are always compatible
  public static final EdgeAttributes COMPATIBLE_EDGE_ATTRIBUTES = new EdgeAttributes() {

    public boolean isCompatible(final EdgeAttributes attr) {

      return equals(attr);
    }
  };
  // a class of node attributes that are always compatible
  public static final NodeAttributes COMPATIBLE_NODE_ATTRIBUTES = new NodeAttributes() {

    public boolean isCompatible(NodeAttributes attr) {

      return equals(attr);
    }
  };
  // default allocation space is roughly 3% of maximum (1/32)
  public static int DYN_ALLOC_NUM = 1;
  public static int DYN_ALLOC_DEN = 32;
  public static boolean DYNAMIC_ALLOCATION = true;
  // default coefficient for a graph
  public static double DEFAULT_GRAPH_COEFFICIENT = 1.0;

  // **************************************************
  //
  // Edges related methods.
  //
  // **************************************************
  private static EdgesMetadata simpleEdgesMetadata(double coefficient) {

    return new SimpleEdgesMetadata(coefficient);
  }

  public static Edges simpleEdges(final EdgesRepresentation rep) {

    return new SimpleEdges(rep, simpleEdgesMetadata(DEFAULT_GRAPH_COEFFICIENT));
  }

  /**
   * Each graph string output by the modified nauty implementation is either an
   * upper triangle or an adjacency matrix representation of the graph. Hence,
   * the client code must make sure the flags passed to nauty are consistent
   * with the flags used to create the representation object.
   */
  public static Edges nautyEdges(final EdgesRepresentation rep,
      double coefficient) {

    return new SimpleEdges(rep, simpleEdgesMetadata(coefficient));
  }


  public static Edges complementEdges(final EdgesRepresentation rep) {

    return new SimpleEdges(rep.complement(), simpleEdgesMetadata(DEFAULT_GRAPH_COEFFICIENT));
  };

  /**
   * Edges attributes.
   */
  public static EdgeAttributes defaultEdgeAttributes() {

    return COMPATIBLE_EDGE_ATTRIBUTES;
  }

  // **************************************************
  //
  // Nodes related methods.
  //
  // **************************************************
  public static Nodes defaultNodes(byte nodeCount) {

    return new SimpleNodes(nodeCount);
  }

  public static Nodes emptyNodes() {

    return new SimpleNodes((byte) 0);
  }

  // **************************************************
  //
  // GraphFamily helper methods.
  //
  // **************************************************
  private static int maxInstances(byte numNodes) {

    if (numNodes == 1) {
      return 0x00000001;
    }
    return 0x00000001 << (numNodes * (numNodes - 1) / 2);
  }

  private static List<Edges> createEdgesList(byte numNodes) {

    int listSize = maxInstances(numNodes);
    if (DYNAMIC_ALLOCATION) {
      listSize = listSize * DYN_ALLOC_NUM / DYN_ALLOC_DEN;
    }
    return new ArrayList<Edges>(listSize);
  }

  // **************************************************
  //
  // GraphFamily instances generated naively.
  //
  // **************************************************
  public static GraphSet naiveGraphSet(final Nodes nodes,
      final EdgesFilter filter) {

    return naiveGraphSet(nodes, filter, false);
  }

  public static GraphSet naiveGraphSet(final Nodes nodes,
      final EdgesFilter filter, boolean isomorphFree) {

    EdgesRepresentationFactory factory = EdgesRepresentationFactory
        .getFactory(nodes.count());
    List<Edges> edgesList = createEdgesList(nodes.count());
    EdgesFilter ef = filter;
    if (isomorphFree) {
      EdgesFilter isof = (new FilterFactory()).isomorphismFilter(nodes,
          edgesList);
      if (ef != null) {
        ef.chain(isof);
      }
      else {
        ef = isof;
      }
    }
    EdgesGenerator generator = new NaiveEdgesGenerator(factory, ef);
    int i = 0;
    Edges e = generator.next(edgesList);
    while (e != null) {
      i++;
      edgesList.add(e);
      e = generator.next(edgesList);
    }
    SimpleGraphSet result = new SimpleGraphSet(nodes, edgesList);
    result.setTags(generator.getTags());
    return result;
  }

  // **************************************************
  //
  // GraphFamily instances generated using nauty.
  //
  // **************************************************
  public static GraphSet nautyGraphSet(final Nodes nodes,
      final EdgesFilter filter, final ProcessWrapper nauty) {

    EdgesRepresentationFactory factory = EdgesRepresentationFactory.getFactory(
        ((NautyInfo) nauty.getProcessInfo()).isUpperTriangle(), nodes.count());
    EdgesGenerator generator = new NautyEdgesGenerator(factory, nauty, filter);
    List<Edges> edgesList = createEdgesList(nodes.count());
    Edges e = generator.next(edgesList);
    while (e != null) {
      edgesList.add(e);
      e = generator.next(edgesList);
    }
    SimpleGraphSet result = new SimpleGraphSet(nodes, edgesList);
    result.setTags(generator.getTags());
    return result;
  }

  /**
   * Because a new constructor is needed for this Nodes class, a nested class
   * had to be defined instead of using an anonymous class, which does not
   * support non-default constructors.
   * 
   * @author Demian Lessa
   */
  static class SimpleNodes implements Nodes {

    byte nodeCount = 0;

    public SimpleNodes(byte nodeCount) {

      this.nodeCount = nodeCount;
    }

    public byte count() {

      return nodeCount;
    }

    public NodeAttributes getAttributes(int nodeID) {

      return COMPATIBLE_NODE_ATTRIBUTES;
    }
  }
}