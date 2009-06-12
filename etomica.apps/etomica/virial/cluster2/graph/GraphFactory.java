package etomica.virial.cluster2.graph;

import java.util.ArrayList;
import java.util.List;

import etomica.virial.cluster2.graph.impl.AbstractEdgesFilter;
import etomica.virial.cluster2.graph.impl.AbstractNodes;
import etomica.virial.cluster2.graph.impl.NaiveEdgesGenerator;
import etomica.virial.cluster2.graph.impl.NautyEdgesGenerator;
import etomica.virial.cluster2.graph.impl.SimpleEdges;
import etomica.virial.cluster2.graph.impl.SimpleEdgesMetadata;
import etomica.virial.cluster2.graph.impl.SimpleGraphSet;
import etomica.virial.cluster2.nauty.NautyInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;

public class GraphFactory {

  // default representation is adjacency matrix
  public static boolean USE_UPPER_TRIANGLE = false;
  // filter flags
  public static final String TAG_RANGE_FILTER = "EdgesRange";
  public static final String TAG_FALSE_FILTER = "FALSE";
  public static final String TAG_TRUE_FILTER = "TRUE";
  public static final String TAG_CONNECTED = "Connected";
  // a class of edge attributes that are always compatible
  public static final EdgeAttributes COMPATIBLE_EDGE_ATTRIBUTES = new EdgeAttributes() {

    @Override
    public boolean isCompatible(final EdgeAttributes attr) {

      return equals(attr);
    }
  };
  // a class of node attributes that are always compatible
  public static final NodeAttributes COMPATIBLE_NODE_ATTRIBUTES = new NodeAttributes() {

    @Override
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

    return new AbstractNodes(nodeCount) {

      @Override
      public NodeAttributes getAttributes(int nodeID) {

        return COMPATIBLE_NODE_ATTRIBUTES;
      }
    };
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

    EdgesRepresentationFactory factory = EdgesRepresentationFactory
        .getFactory(nodes.count());
    EdgesGenerator generator = new NaiveEdgesGenerator(factory, filter);
    List<Edges> edgesList = createEdgesList(nodes.count());
    while (generator.hasNext()) {
      edgesList.add(generator.next());
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
    while (generator.hasNext()) {
      edgesList.add(generator.next());
    }
    SimpleGraphSet result = new SimpleGraphSet(nodes, edgesList);
    result.setTags(generator.getTags());
    return result;
  }

  // **************************************************
  //
  // EdgesFilter instances.
  //
  // **************************************************
  public static EdgesFilter trueFilter() {

    return new AbstractEdgesFilter() {

      @Override
      protected boolean doAccept(Edges edges) {

        return true;
      }

      @Override
      protected String tag() {

        return TAG_TRUE_FILTER;
      }
    };
  }

  public static EdgesFilter falseFilter() {

    return new AbstractEdgesFilter() {

      @Override
      protected boolean doAccept(Edges edges) {

        return false;
      }

      @Override
      protected String tag() {

        return TAG_FALSE_FILTER;
      }
    };
  }

  public static EdgesFilter connectedFilter(final Nodes nodes) {

    return new AbstractEdgesFilter() {

      @Override
      protected boolean doAccept(Edges edges) {

        return Algorithms.isConnected(nodes, edges);
      }

      @Override
      protected String tag() {

        return TAG_CONNECTED;
      }
    };
  }

  public static EdgesFilter rangeFilter(int minEdges, int maxEdges) {

    return new RangeFilter(minEdges, maxEdges);
  }

  static class RangeFilter extends AbstractEdgesFilter {

    private static int maxEdges;
    private static int minEdges;

    RangeFilter(int min, int max) {

      minEdges = min;
      maxEdges = max;
    }

    @Override
    protected boolean doAccept(Edges edges) {

      int count = edges.count();
      return (count >= minEdges) && (count <= maxEdges);
    }

    @Override
    protected String tag() {

      return TAG_RANGE_FILTER + ":" + minEdges + ":" + maxEdges;
    }
  }
}