package etomica.graph.model.impl;

import java.util.ArrayList;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import etomica.graph.model.Bitmap;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Coefficient;
import etomica.graph.model.Edge;
import etomica.graph.model.EdgeVisitor;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Node;
import etomica.graph.model.NodeVisitor;
import etomica.graph.traversal.DepthFirst;
import etomica.graph.traversal.Traversal;

import static etomica.graph.model.Metadata.*;

public class GraphImpl implements Graph {

  private Bitmap store;
  private Node[] nodes;
  private Coefficient coefficient;
  private Map<Byte, Edge> edges = new HashMap<Byte, Edge>();

  private GraphImpl(Node[] nodes, Map<Byte, Edge> edges, Bitmap store, Coefficient coefficient) {

    this.nodes = nodes;
    this.edges = edges;
    this.store = store;
    this.coefficient = coefficient;
  }

  public GraphImpl(Node[] nodes) {

    this.coefficient = GraphFactory.createCoefficient();
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
    this.nodes = new Node[nodeCount];
    for (byte i = 0; i < this.nodes.length; i++) {
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
        // check store
        int storeOrder = store.compareTo(other.getStore());
        if (storeOrder != 0) {
          return storeOrder;
        }
        // check nodes
        for (int nodeId = 0; nodeId < nodes().size(); nodeId++) {
          int nodeOrder = nodes().get(nodeId).compareTo(other.nodes().get(nodeId));
          if (nodeOrder != 0) {
            return nodeOrder;
          }
        }
        // check edges
        for (byte edgeId = 0; edgeId < edges().size(); edgeId++) {
          Edge edge = edges.get(edgeId);
          if (edge != null) {
            int edgeOrder = edge.compareTo(other.getEdge(edge.getId()));
            if (edgeOrder != 0) {
              return edgeOrder;
            }
          }
        }
        return 0;
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

    return new ArrayList<Edge>(edges.values());
  }

  public String edgesToString() {

    // use curly brackets for a set
    String result = "";
    // boolean isMono = (getColors(TYPE_EDGE_ANY).size() == 1);
    boolean first = true;
    for (Byte edgeId : edges.keySet()) {
      if (!first) {
        result += ", ";
      }
      result += "(" + getFromNode(edgeId) + "," + getToNode(edgeId) + ")";
      // if (!isMono) {
      result += edges.get(edgeId).getColor();
      // }
      first = false;
    }
    return "{" + result + "}";
  }

  @Override
  public boolean equals(Object obj) {

    if (obj instanceof Graph) {
      Graph other = (Graph) obj;
      // check coefficient
      if (!coefficient.equals(other.coefficient())) {
        return false;
      }
      // check node count
      if (nodes.length != other.nodeCount()) {
        return false;
      }
      // check edge count
      if (edges.size() != other.edgeCount()) {
        return false;
      }
      // check nodes
      for (int nodeId = 0; nodeId < nodes.length; nodeId++) {
        if (!nodes[nodeId].equals(other.nodes().get(nodeId))) {
          return false;
        }
      }
      // check edges
      for (Edge edge : edges.values()) {
        if (!edge.equals(other.getEdge(edge.getId()))) {
          return false;
        }
      }
      // no need to check the store because the edgeId values
      // implicitly test for the bits set in the store bitmap
      return true;
    }
    return false;
  }

  // computes the signature of the graph:
  // - add a string CCStr for the number of connected components
  // - add a string RStr for each root node, in ID order
  // - add a string FStr for each field node color, in color order
  // - add a string EStr for each edge node color, in color order
  // - CCStr = '/CC' || count of connected components
  // - RStr = '/R' || (color || Id || outDegree)* ordered by Id
  // - FStr = '/F' || (color count)* ordered by color || (outDegree)* ordered
  // - EStr = '/E' || (color count)* ordered by color
  public String getSignature() {
    String result = "/R";
    SortedMap<Character, List<Byte>> fieldColorMap = new TreeMap<Character, List<Byte>>();
    byte rootCount = 0;
    for (byte nodeId = 0; nodeId < nodes.length; nodeId++) {
      if (nodes[nodeId].getType() == TYPE_NODE_ROOT) {
        rootCount++;
      }
    }
    result += "" + rootCount; 
    for (byte nodeId = 0; nodeId < nodes.length; nodeId++) {
      List<Byte> degrees = fieldColorMap.get(nodes[nodeId].getColor());
      if (degrees == null) {
        degrees = new ArrayList<Byte>();
        fieldColorMap.put(nodes[nodeId].getColor(), degrees);
      }
      degrees.add(getOutDegree(nodeId));
    }
    result += "/F";
    for (Character color : fieldColorMap.keySet()) {
      List<Byte> degrees = fieldColorMap.get(color);
      Collections.sort(degrees);
      result += "" + color;
      for (Byte degree : degrees) {
        result += ":" + degree;
      }
    }
    SortedMap<Character, Byte> edgeColorMap = new TreeMap<Character, Byte>();
    for (Edge edge : edges.values()) {
      Byte count = edgeColorMap.get(edge.getColor());
      if (count == null) {
        edgeColorMap.put(edge.getColor(), (byte) 1);
      }
      else {
        edgeColorMap.put(edge.getColor(), (byte) (count + 1));
      }
    }
    result += "/E";
    for (Character color : edgeColorMap.keySet()) {
      result += "" + color + edgeColorMap.get(color);
    }
    Traversal t = new DepthFirst();
    byte cc = t.traverseAll(this, null);
    result = "/CC" + cc + result;
    return result;
  }

  private Set<Character> getColors(char type) {

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

  public byte edgeCount() {

    return (byte) store.bitCount();
  }

  // (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),(0,4),(1,4),(2,4),(3,4)...,(n-1,n)
  // edgeId = (toNode)*(toNode-1)/2 + fromNode
  public byte getEdgeId(byte fromNode, byte toNode) {

    assert (fromNode != toNode);
    if (fromNode > toNode) {
      return getEdgeId(toNode, fromNode);
    }
    return (byte) (fromNode + toNode * (toNode - 1) / 2);
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
  }

  public Node getNode(byte node) {

    return nodes[node];
  }

  public byte nodeCount() {

    return (byte) this.nodes.length;
  }

  public byte getOutDegree(byte node) {

    byte result = 0;
    for (byte i = 0; i < nodeCount(); i++) {
      if (i != node) {
        result += hasEdge(node, i) ? 1 : 0;
      }
    }
    return result;
  }

  public byte getOutNode(byte node, byte index) {

    byte found = 0;
    for (byte i = 0; i < nodeCount(); i++) {
      if (i != node) {
        found += hasEdge(node, i) ? 1 : 0;
      }
      if (index == found - 1) {
        return i;
      }
    }
    return (byte) 0xFF;
  }

  public Bitmap getStore() {

    return store;
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
  }

  public boolean hasEdge(byte fromNode, byte toNode) {

    return store.testBit(getEdgeId(fromNode, toNode));
  }

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

  public Edge putEdge(byte fromNode, byte toNode) {

    byte edgeId = getEdgeId(fromNode, toNode);
    store.setBit(edgeId);
    Edge edge = GraphFactory.createEdge(edgeId);
    edges.put(edgeId, edge);
    return edge;
  }

  public void putEdge(byte edgeId) {

    store.setBit(edgeId);
    edges.put(edgeId, GraphFactory.createEdge(edgeId));
  }

  @Override
  public String toString() {

    return coefficient.toString() + /* " :: " + getSignature() + */" :: " + nodesToString() + " :: "
        + edgesToString();
  }

  public String toSVG(int dim) {

    int nodeR = dim / 12;
    int oX = dim / 2;
    int oY = dim / 2;
    int graphR = dim / 2 - nodeR;

    // first try 0, then A, then use this fixed one-- should work for 0
    double oA = Math.PI / 2;
    double angle = 2 * Math.PI / nodes.length;

    double[] x = new double[nodes.length];
    double[] y = new double[nodes.length];

    String svgNodes = "";
    for (int i = 0; i < nodes.length; i++) {
      x[i] = oX + graphR * Math.cos(oA + i * angle);
      y[i] = oY + graphR * Math.sin(oA + i * angle);
      if (nodes[i].getType() == TYPE_NODE_ROOT) {
        svgNodes += String.format(
            "<circle style=\"fill: white; stroke: %s\" r=\"%d\" cx=\"%.2f\" cy=\"%.2f\"/>\n",
            COLORS[COLOR_CODES.indexOf(nodes[i].getColor())], nodeR, x[i], y[i]);
      }
      else {
        svgNodes += String.format(
            "<circle style=\"fill: %s; stroke: black\" r=\"%d\" cx=\"%.2f\" cy=\"%.2f\"/>\n",
            COLORS[COLOR_CODES.indexOf(nodes[i].getColor())], nodeR, x[i], y[i]);
      }
    }
    String svgEdges = "";
    for (Edge e : edges()) {
      byte nodeFrom = getFromNode(e.getId());
      byte nodeTo = getToNode(e.getId());
      svgEdges += String.format(
          "<line style=\"stroke: %s\" x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>\n",
          COLORS[COLOR_CODES.indexOf(e.getColor())], x[nodeFrom], y[nodeFrom], x[nodeTo], y[nodeTo]);
    }
    return svgEdges + svgNodes;
  }

  public void visitEdges(EdgeVisitor visitor) {

    for (Edge edge : edges.values()) {
      if (!visitor.visit(edge)) {
        return;
      }
    }
  }

  public void visitNodes(NodeVisitor visitor) {

    for (byte i = 0; i < nodeCount(); i++) {
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
