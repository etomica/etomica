package etomica.virial.cluster2.graph;

import java.util.ArrayList;
import java.util.Collections;
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
  // all field nodes have the same attributes and are compatible amongst each
  // other
  public static final NodeAttributes FIELD_NODE_ATTRIBUTES = new NodeAttributes() {

    public boolean isCompatible(NodeAttributes attr) {

      return isSameColor(attr);
    }

    public boolean isSameColor(NodeAttributes attr) {

      return equals(attr);
    }
  };
  // template root node attributes
  public static final NodeAttributes ROOT_NODE_ATTRIBUTES = new RootNodeAttributes(
      (byte) 0);
  // default allocation space is roughly 3% of maximum (1/32)
  public static int DYN_ALLOC_NUM = 1;
  public static int DYN_ALLOC_DEN = 32;
  public static boolean DYNAMIC_ALLOCATION = true;
  // default coefficient for a graph
  public static double DEFAULT_GRAPH_COEFFICIENT = 1.0;

  // **************************************************
  //
  // A graph is simply a container for nodes and edges,
  // and a GraphSet is a container for nodes and a set
  // of edges.
  //
  // **************************************************
  public static Graph simpleGraph(final Nodes nodes, final Edges edges) {

    return new Graph() {

      public Edges getEdges() {

        return edges;
      }

      public Nodes getNodes() {

        return nodes;
      }
    };
  }

  public static GraphSet simpleGraphSet(final Nodes nodes,
      final List<Edges> list) {

    return new SimpleGraphSet(nodes, list);
  }

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

  public static Edges complementEdges(final Edges edges) {

    return new SimpleEdges(edges.getRepresentation().complement(),
        simpleEdgesMetadata(edges.getMetadata().getCoefficient()));
  };

  public static Edges canonicalEdges(final Edges edges) {

    return new SimpleEdges(edges.getRepresentation().canonical(),
        simpleEdgesMetadata(edges.getMetadata().getCoefficient()));
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

    return new FieldNodes(nodeCount);
  }

  public static Nodes defaultNodes(byte fieldNodes, byte rootNodes) {

    return new SimpleNodes(fieldNodes, rootNodes);
  }

  public static Nodes emptyNodes() {

    return new FieldNodes((byte) 0);
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

  public static List<Edges> createEdgesList(byte numNodes) {

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
    SimpleGraphSet result = (SimpleGraphSet)simpleGraphSet(nodes, edgesList);
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
    SimpleGraphSet result = (SimpleGraphSet)simpleGraphSet(nodes, edgesList);
    result.setTags(generator.getTags());
    return result;
  }

  /**
   * Because a new constructor is needed for this Nodes class, a nested class
   * had to be defined instead of using an anonymous class, which does not
   * support non-default constructors. This class supports only field nodes, all
   * of which are compatible with each other.
   * 
   * @author Demian Lessa
   */
  static class FieldNodes implements Nodes {

    private byte nodeCount = 0;
    private List<Integer> partition;

    public FieldNodes(byte nodeCount) {

      assert (nodeCount > 0);
      this.nodeCount = nodeCount;
      List<Integer> p = new ArrayList<Integer>(nodeCount);
      for (int i = 0; i < nodeCount; i++) {
        p.add(i);
      }
      partition = Collections.unmodifiableList(p);
    }

    public byte count() {

      return nodeCount;
    }

    /**
     * All field nodes are compatible with each other.
     */
    public NodeAttributes getAttributes(int nodeID) {

      return FIELD_NODE_ATTRIBUTES;
    }

    public List<Integer> getPartition(int partitionID) {

      assert (partitionID == 0);
      return partition;
    }

    public int getPartitionCount() {

      return 1;
    }
  }

  /**
   * Because a new constructor is needed for this Nodes class, a nested class
   * had to be defined instead of using an anonymous class, which does not
   * support non-default constructors. This class supports both field and root
   * nodes. Field nodes are compatible with field nodes, and root nodes are
   * incompatible with all nodes. We assume that the first fieldNodeCount nodes
   * are field nodes, and the remaining are root nodes.
   * 
   * @author Demian Lessa
   */
  static class SimpleNodes implements Nodes {

    private byte fieldNodeCount = 0;
    private byte rootNodeCount = 0;
    private NodeAttributes[] rootNodeAttributes;
    private List<List<Integer>> partition;

    public SimpleNodes(byte fieldNodeCount, byte rootNodeCount) {

      assert (fieldNodeCount >= 0);
      assert (rootNodeCount >= 0);
      assert (fieldNodeCount + rootNodeCount > 0);
      assert (fieldNodeCount + rootNodeCount <= 255);
      this.fieldNodeCount = fieldNodeCount;
      this.rootNodeCount = rootNodeCount;
      this.rootNodeAttributes = new RootNodeAttributes[rootNodeCount];
      for (byte i = 0; i < rootNodeCount; i++) {
        rootNodeAttributes[i] = new RootNodeAttributes(i);
      }
      int partitionCount = rootNodeCount + (fieldNodeCount > 0 ? 1 : 0);
      List<List<Integer>> ps = new ArrayList<List<Integer>>(partitionCount);
      for (int i = 0; i < rootNodeCount; i++) {
        ArrayList<Integer> p = new ArrayList<Integer>(1);
        p.add(i);
        ps.add(Collections.unmodifiableList(p));
      }
      if (fieldNodeCount > 0) {
        ArrayList<Integer> p = new ArrayList<Integer>(fieldNodeCount);
        for (int i = 0; i < fieldNodeCount; i++) {
          p.add(rootNodeCount + i);
        }
        ps.add(Collections.unmodifiableList(p));
      }
      partition = Collections.unmodifiableList(ps);
    }

    public byte count() {

      return (byte) (fieldNodeCount + rootNodeCount);
    }

    // from 0..fieldNodeCount-1 => field node
    // from fieldNodeCount..count() => root node
    public NodeAttributes getAttributes(int nodeID) {

      assert (nodeID >= 0 && nodeID < count());
      if (nodeID < rootNodeCount) {
        return rootNodeAttributes[nodeID];
      }
      else {
        return FIELD_NODE_ATTRIBUTES;
      }
    }

    public List<Integer> getPartition(int partitionID) {

      assert (partitionID >= 0 && partitionID < getPartitionCount());
      return partition.get(partitionID);
    }

    public int getPartitionCount() {

      return partition.size();
    }
  }

  static class RootNodeAttributes implements NodeAttributes {

    private byte nodeID;

    public RootNodeAttributes(byte nodeID) {

      this.nodeID = nodeID;
    }

    public boolean isCompatible(NodeAttributes attr) {

      return isSameColor(attr) && ((RootNodeAttributes) attr).nodeID == nodeID;
    }

    public boolean isSameColor(NodeAttributes attr) {

      return (attr instanceof RootNodeAttributes);
    }
  }
}