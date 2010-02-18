package etomica.graph.model.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.graph.model.Bitmap;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Coefficient;
import etomica.graph.model.Edge;
import etomica.graph.model.EdgeVisitor;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Node;
import etomica.graph.model.NodeVisitor;

import static etomica.graph.model.Metadata.*;

public class GraphImpl implements Graph {

  private Bitmap store;
  private byte[] labels;
  private Node[] nodes;
  private Coefficient coefficient;
  private Map<Byte, Edge> edges = new HashMap<Byte, Edge>();

  public GraphImpl(Node[] nodes, Map<Byte, Edge> edges, Bitmap store, Coefficient coefficient) {

    this.nodes = nodes;
    this.edges = edges;
    this.store = store;
    this.coefficient = coefficient;
    this.labels = new byte[nodes.length];
    for (byte i = 0; i < nodes.length; i++) {
      this.labels[i] = i;
    }
  }

  public GraphImpl(Node[] nodes) {

    this.coefficient = GraphFactory.createCoefficient();
    this.labels = new byte[nodes.length];
    this.nodes = new Node[nodes.length];
    for (byte i = 0; i < nodes.length; i++) {
      this.labels[i] = i;
    }
    this.nodes = nodes;
    this.store = BitmapFactory.createBitmap((byte) nodes.length, false);
    createEdges();
  }

  public GraphImpl(byte nodeCount) {

    this(nodeCount, (byte) 0, BitmapFactory.createBitmap(nodeCount, false));
  }

  public GraphImpl(byte nodeCount, byte rootNodeCount) {

    this(nodeCount, rootNodeCount, BitmapFactory.createBitmap(nodeCount, false));
  }

  public GraphImpl(byte nodeCount, Bitmap store) {

    this(nodeCount, (byte) 0, store, GraphFactory.createCoefficient());
  }

  public GraphImpl(byte nodeCount, byte rootNodeCount, Bitmap store) {

    this(nodeCount, rootNodeCount, store, GraphFactory.createCoefficient());
  }

  public GraphImpl(byte nodeCount, byte rootNodeCount, Bitmap store, Coefficient coefficient) {

    this.coefficient = coefficient;
    this.labels = new byte[nodeCount];
    this.nodes = new Node[nodeCount];
    for (byte i = 0; i < this.nodes.length; i++) {
      this.labels[i] = i;
      this.nodes[i] = GraphFactory.createNode(i, i < rootNodeCount);
    }
    this.store = store;
    createEdges();
  }

  public Coefficient coefficient() {

    return this.coefficient;
  }

  public int compareTo(Graph other) {

    // two graphs are indistinguishable if and only if they have the same number of nodes,
    // the same number of edges, the same edges bitmap representation (store), and the
    // same node color strings; the ordering based on color strings is arbitrary, yet
    // useful as a tie breaker; the actual colors mapped to each of these abstract
    // colors can be freely defined
    if (nodes().size() == other.nodes().size()) {
      if (edges.size() == other.edges().size()) {
        int result = store.compareTo(other.getStore());
        if (result == 0) {
          for (int nodeId = 0; nodeId < nodes().size(); nodeId++) {
            if (nodes().get(nodeId).getColor() < other.nodes().get(nodeId).getColor()) {
              return -1;
            }
            else if (nodes().get(nodeId).getColor() > other.nodes().get(nodeId).getColor()) {
              return 1;
            }
          }
        }
        return result;
      }
      else if (edges().size() < other.edges().size()) {
        return -1;
      }
      else {
        return 1;
      }
    }
    else if (nodes().size() < other.nodes().size()) {
      return -1;
    }
    else {
      return 1;
    }
  }

  public Graph copy() {

    // copy the nodes
    Node[] nodesCopy = new Node[nodes.length];
    for (int nodeId = 0; nodeId < nodes.length; nodeId++) {
      nodesCopy[nodeId] = nodes[nodeId].copy();
    }
    // copy the edges
    Map<Byte, Edge> edgesCopy = new HashMap<Byte, Edge>();
    for (Byte edgeId : edges.keySet()) {
      edgesCopy.put(edgeId, edges.get(edgeId).copy());
    }
    return new GraphImpl(nodesCopy, edgesCopy, store.copy(), coefficient.copy());
  }

  protected void createEdges() {

    for (byte i = 0; i < nodes.length; i++) {
      for (byte j = (byte) (i + 1); j < nodes.length; j++) {
        if (hasEdge(i, j)) {
          byte edgeId = getEdgeId(i, j);
          edges.put(edgeId, GraphFactory.createEdge(edgeId));
        }
      }
    }
  }

  public void deleteEdge(byte fromNode, byte toNode) {

    byte edgeId = getEdgeId(fromNode, toNode);
    edges.remove(edgeId);
    store.clearBit(edgeId);
  }

  public List<Edge> edges() {

    EdgeCollectorVisitor v = new EdgeCollectorVisitor();
    visitEdges(v);
    return v.getEdges();
  }

  public String edgesToString() {

    // use curly brackets for a set
    String result = "";
    boolean isMono = (getColors(TYPE_EDGE_ANY).size() == 1);
    boolean first = true;
    for (Byte edgeId : edges.keySet()) {
      if (!first) {
        result += ", ";
      }
      result += "(" + getFromNode(edgeId) + "," + getToNode(edgeId) + ")";
      if (!isMono) {
        result += edges.get(edgeId).getColor();
      }
      first = false;
    }
    return "{" + result + "}";
  }

  @Override
  public boolean equals(Object obj) {

    if (obj instanceof Graph) {
      Graph other = (Graph) obj;
      return toString().equals(other.toString());
    }
    return false;
  }

  public Set<Character> getColors(char type) {

    if (type == TYPE_EDGE_ANY) {
      EdgeColorVisitor v = new EdgeColorVisitor();
      visitEdges(v);
      return v.getColors();
    }
    else {
      NodeColorVisitor v = new NodeColorVisitor(type);
      visitNodes(v);
      return v.getColors();
    }
  }

  public Edge getEdge(byte fromNode, byte toNode) {

    return edges.get(getEdgeId(fromNode, toNode));
  }

  public Edge getEdge(byte edgeId) {

    return edges.get(edgeId);
  }

  public byte getEdgeCount() {

    return (byte) store.bitCount();
  }

  // (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),(0,4),(1,4),(2,4),(3,4)...,(n-1,n)
  // edgeId = (toNode)*(toNode-1)/2 + fromNode
  // (0,1) => (1)(0)/2 + 0 = 0
  // (0,2) => (2)(1)/2 + 0 = 1
  // (1,2) => (2)(1)/2 + 1 = 2
  // (0,3) => (3)(2)/2 + 0 = 3
  // (1,3) => (3)(2)/2 + 1 = 4
  // (2,3) => (3)(2)/2 + 2 = 5
  // (0,4) => (4)(3)/2 + 0 = 6

  public byte getEdgeId(byte fromNode, byte toNode) {

    assert (fromNode != toNode);
    if (fromNode > toNode) {
      return getEdgeId(toNode, fromNode);
    }
    return (byte) (fromNode + toNode * (toNode - 1) / 2);
    // return (byte) ((toNode - fromNode - 1) + sumMaxEdges((byte) (fromNode - 1)));
  }

  public byte getFromLabel(byte edge) {

    return getLabel(getFromNode(edge));
  }

  public byte getFromNode(byte edge) {

    for (int toNode = 1; toNode < nodes.length; toNode++) {
      // 0, 1, 3, 6, 10, 15, ...
      byte sectionStart = (byte) (toNode * (toNode - 1) / 2);
      // edge - candidate is the offset of the fromNode in a toNode section
      if (edge - sectionStart < toNode) {
        return (byte) (edge - sectionStart);
      }
    }
    return 0;
//    byte offset = (byte) (maxEdges(fromNode) - 1);
//    while (edge > offset) {
//      fromNode++;
//      offset += maxEdges(fromNode);
//    }
//    return fromNode;
  }

  // returns the label of the given node
  public byte getLabel(byte node) {

    for (byte i = 0; i < labels.length; i++) {
      if (labels[i] == node) {
        return i;
      }
    }
    assert (false);
    return (byte) 0xFF;
  }

  public byte getLabeledNode(byte label) {

    return labels[label];
  }

  public Node getNode(byte node) {

    return nodes[node];
  }

  public byte getNodeCount() {

    return (byte) this.nodes.length;
  }

  public byte getOutDegree(byte node) {

    byte result = 0;
    for (byte i = 0; i < getNodeCount(); i++) {
      if (i != node) {
        result += hasEdge(node, i) ? 1 : 0;
      }
    }
    return result;
  }

  public byte getOutNode(byte node, byte index) {

    byte found = 0;
    for (byte i = 0; i < getNodeCount(); i++) {
      if (i != node) {
        found += hasEdge(node, i) ? 1 : 0;
      }
      if (index == found - 1) {
        return i;
      }
    }
    return (byte) 0xFF;
  }

  public List<Byte> getPartition(char type, char color) {

    if (type == TYPE_EDGE_ANY) {
      EdgeColorPartitionVisitor v = new EdgeColorPartitionVisitor(color);
      visitEdges(v);
      return v.getPartition();
    }
    else {
      NodeColorPartitionVisitor v = new NodeColorPartitionVisitor(type, color);
      visitNodes(v);
      return v.getPartition();
    }
  }

  public Bitmap getStore() {

    return store;
  }

  public byte getToLabel(byte edge) {

    return getLabel(getToNode(edge));
  }

  public byte getToNode(byte edge) {

    for (int toNode = 1; toNode < nodes.length; toNode++) {
      // 0, 1, 3, 6, 10, 15, ...
      byte sectionStart = (byte) (toNode * (toNode - 1) / 2);
      // edge - candidate is the offset of the fromNode in a toNode section
      if (edge - sectionStart < toNode) {
        return (byte) toNode;
      }
    }
    return 0;
//    byte fromNode = getFromNode(edge);
//    int fromEdge = getEdgeId(fromNode, (byte) (fromNode + 1));
//    return (byte) ((fromNode + 1) + (edge - fromEdge));
  }

  public boolean hasEdge(byte fromNode, byte toNode) {

    return store.testBit(getEdgeId(fromNode, toNode));
  }

//  protected byte maxEdges(byte node) {
//
//    return (byte) (getNodeCount() - node - 1);
//  }

  public List<Node> nodes() {

    return Arrays.asList(nodes);
  }

  public String nodesToString() {

    // use square brackets for a list
    String result = "[";
    for (byte nodeId = 0; nodeId < nodes.length; nodeId++) {
      result += nodes[nodeId];
      if (nodeId < nodes.length - 1) {
        result += ", ";
      }
    }
    result += "]";
    return result;
  }

  public void putEdge(byte fromNode, byte toNode) {

    byte edgeId = getEdgeId(fromNode, toNode);
    store.setBit(edgeId);
    edges.put(edgeId, GraphFactory.createEdge(edgeId));
  }

  public void putEdge(byte edgeId) {

    store.setBit(edgeId);
    edges.put(edgeId, GraphFactory.createEdge(edgeId));
  }

  // use this label to refer to this node
  public void setLabel(byte label, byte node) {

    labels[label] = node;
  }

//  private byte sumMaxEdges(byte lastNode) {
//
//    byte result = 0;
//    for (byte node = 0; node <= lastNode; node++) {
//      result += maxEdges(node);
//    }
//    return result;
//  }

  @Override
  public String toString() {

    return coefficient.toString() + " :: " + nodesToString() + " :: " + edgesToString();
  }

  public void visitEdges(EdgeVisitor visitor) {

    for (Edge edge : edges.values()) {
      if (!visitor.visit(edge)) {
        return;
      }
    }
  }

  public void visitNodes(NodeVisitor visitor) {

    for (byte i = 0; i < getNodeCount(); i++) {
      if (!visitor.visit(getNode(i))) {
        return;
      }
    }
  }
}

class EdgeCollectorVisitor implements EdgeVisitor {

  List<Edge> edges = new ArrayList<Edge>();

  public List<Edge> getEdges() {

    return edges;
  }

  public boolean visit(Edge edge) {

    edges.add(edge);
    return true;
  }
}

class EdgeColorPartitionVisitor implements EdgeVisitor {

  List<Byte> partition = new ArrayList<Byte>();
  private char color;

  public EdgeColorPartitionVisitor(char color) {

    this.color = color;
  }

  public List<Byte> getPartition() {

    return partition;
  }

  public boolean visit(Edge edge) {

    if (edge.getColor() == color) {
      partition.add(edge.getId());
    }
    return true;
  }
}

class EdgeColorVisitor implements EdgeVisitor {

  Set<Character> colors = new HashSet<Character>();

  public Set<Character> getColors() {

    return colors;
  }

  public boolean visit(Edge edge) {

    colors.add(edge.getColor());
    return true;
  }
}

class NodeColorPartitionVisitor implements NodeVisitor {

  List<Byte> partition = new ArrayList<Byte>();
  private char color;
  private char type;

  public NodeColorPartitionVisitor(char type, char color) {

    this.type = type;
    this.color = color;
  }

  public List<Byte> getPartition() {

    return partition;
  }

  public boolean visit(Node node) {

    if (node.getType() == type && node.getColor() == color) {
      partition.add(node.getId());
    }
    return true;
  }
}

class NodeColorVisitor implements NodeVisitor {

  Set<Character> colors = new HashSet<Character>();
  private char type;

  public NodeColorVisitor(char type) {

    this.type = type;
  }

  public Set<Character> getColors() {

    return colors;
  }

  public boolean visit(Node node) {

    if (node.getType() == type) {
      colors.add(node.getColor());
    }
    return true;
  }
}