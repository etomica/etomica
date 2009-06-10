package etomica.virial.cluster2.graph;

import java.util.ArrayList;
import java.util.List;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.impl.AbstractEdgesFilter;
import etomica.virial.cluster2.graph.impl.AbstractNodes;
import etomica.virial.cluster2.graph.impl.AdjacencyMatrixRepresentation;
import etomica.virial.cluster2.graph.impl.NaiveEdgesGenerator;
import etomica.virial.cluster2.graph.impl.NautyEdgesGenerator;
import etomica.virial.cluster2.graph.impl.SimpleEdges;
import etomica.virial.cluster2.graph.impl.SimpleEdgesMetadata;
import etomica.virial.cluster2.graph.impl.SimpleGraphSet;
import etomica.virial.cluster2.graph.impl.UpperTriangleRepresentation;
import etomica.virial.cluster2.nauty.NautyInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;

public class GraphFactory {

  // filter flags
  public static final String FLAG_NULL = "NULL";
  public static final String FLAG_CONNECTED = "Connected";
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
  // assume graphs are represented by adjacency matrices
  public static boolean useUpperTriangle = false;
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

  public static Edges simpleEdges(final Bitmap edges,
      final EdgesRepresentation decoder) {

    return new SimpleEdges(edges, decoder,
        simpleEdgesMetadata(DEFAULT_GRAPH_COEFFICIENT));
  }

  /**
   * Each graph string output by the modified nauty implementation is either an
   * upper triangle or an adjacency matrix representation of the graph. Hence,
   * the client code must make sure the flags passed to nauty are consistent
   * with the flags used to create the representation object.
   */
  public static Edges nautyEdges(final String edges,
      final EdgesRepresentation rep, double coefficient) {

    return new SimpleEdges(BitmapFactory.getBitmap(edges), rep,
        simpleEdgesMetadata(coefficient));
  }

  /**
   * Edges representation.
   */
  private static EdgesRepresentation getRepresentation(boolean isUpperTriangle,
      byte nodeCount) {

    if (isUpperTriangle) {
      return new UpperTriangleRepresentation(nodeCount);
    }
    else {
      return new AdjacencyMatrixRepresentation(nodeCount);
    }
  }

  private static EdgesRepresentation defaultRepresentation(byte nodeCount) {

    return getRepresentation(useUpperTriangle, nodeCount);
  }

  /**
   * Edges attributes.
   */
  public static EdgeAttributes defaultEdgeAttributes() {

    return COMPATIBLE_EDGE_ATTRIBUTES;
  }

  /**
   * Edges storage.
   */
  private static Bitmap getBitmap(final EdgesRepresentation rep,
      byte nodeCount, boolean isSet) {

    if (nodeCount == 0) {
      return BitmapFactory.EMPTY;
    }
    if (nodeCount == 1) {
      return BitmapFactory.ZERO;
    }
    else {
      int capacity;
      if (rep instanceof AdjacencyMatrixRepresentation) {
        capacity = nodeCount * nodeCount;
      }
      else if (rep instanceof UpperTriangleRepresentation) {
        capacity = nodeCount * (nodeCount - 1) / 2;
      }
      else {
        return null;
      }
      return BitmapFactory.getBitmap(capacity, isSet);
    }
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

    EdgesRepresentation rep = GraphFactory.defaultRepresentation(nodes.count());
    Bitmap first = GraphFactory.getBitmap(rep, nodes.count(), false);
    Bitmap max = GraphFactory.getBitmap(rep, nodes.count(), true);
    EdgesGenerator generator = new NaiveEdgesGenerator(rep, first, max, filter);
    List<Edges> edgesList = createEdgesList(nodes.count());
    while (generator.hasNext()) {
      edgesList.add(generator.next());
    }
    return new SimpleGraphSet(nodes, edgesList);
  }

  // **************************************************
  //
  // GraphFamily instances generated using nauty.
  //
  // **************************************************
  public static GraphSet nautyGraphSet(final Nodes nodes,
      final EdgesFilter filter, final ProcessWrapper nauty) {

    EdgesRepresentation rep = GraphFactory.getRepresentation(nauty
        .getProcessInfo().getTags().contains(NautyInfo.TAG_UPPER_TRIANGLE),
        nodes.count());
    EdgesGenerator generator = new NautyEdgesGenerator(rep, nauty, filter);
    List<Edges> edgesList = createEdgesList(nodes.count());
    while (generator.hasNext()) {
      edgesList.add(generator.next());
    }
    return new SimpleGraphSet(nodes, edgesList);
  }

  // **************************************************
  //
  // EdgesFilter instances.
  //
  // **************************************************
  public static EdgesFilter nullFilter(final EdgesFilter filter) {

    return new AbstractEdgesFilter(filter) {

      @Override
      protected boolean doAccept(Edges edges) {

        return false;
      }

      @Override
      protected String tag() {

        return FLAG_NULL;
      }
    };
  }

  public static EdgesFilter connectedFilter(final Nodes nodes,
      final EdgesFilter filter) {

    return new AbstractEdgesFilter(filter) {

      @Override
      protected boolean doAccept(Edges edges) {

        return Algorithms.isConnected(nodes, edges);
      }

      @Override
      protected String tag() {

        return FLAG_CONNECTED;
      }
    };
  }
}