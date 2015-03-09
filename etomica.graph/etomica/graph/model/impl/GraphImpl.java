/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.impl;

import static etomica.graph.model.Metadata.COLORS;
import static etomica.graph.model.Metadata.COLOR_CODES;
import static etomica.graph.model.Metadata.COLOR_MAP;
import static etomica.graph.model.Metadata.DASH_MAP;
import static etomica.graph.model.Metadata.BOND_ORDER_MAP;
import static etomica.graph.model.Metadata.TYPE_EDGE_ANY;
import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
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

public class GraphImpl implements Graph {
	
  public static boolean useReverseEdges = false;

  private final Bitmap store;
  private final Node[] nodes;
  private List<Node> nodeList;
  private List<Edge> edgeList;
  private final Coefficient coefficient;
  private final Edge[] edges, reverseEdges;
  private int[] factors = new int[0];

  private GraphImpl(Node[] nodes, Edge[] edges, Bitmap store, Coefficient coefficient) {

    this.nodes = nodes;
    this.edges = edges;
    this.store = store;
    this.coefficient = coefficient;
    if (useReverseEdges){
	    reverseEdges = new Edge[edges.length];//for multiple site models
	    createReverseEdges();
    } 
    else {
    	reverseEdges = null;
    }
  }

  public GraphImpl(Node[] nodes) {

    this.coefficient = GraphFactory.createCoefficient();
    this.nodes = nodes;
    this.store = BitmapFactory.createBitmap((byte) nodes.length, false);
    edges = new Edge[nodes.length*(nodes.length-1)/2];
    if (useReverseEdges){
    	reverseEdges = new Edge[edges.length];//for multiple site models
    } 
    else {
    	reverseEdges = null;
    }
    createEdges();
  }
  
  public void createReverseEdges(){//for Wertheim multiple site model
	  if(!useReverseEdges) return;
	  for (byte i = 0; i< edges.length; i++){
		  createReverseEdge(i);
	  }
  }
  
  protected void createReverseEdge (byte i){
	  reverseEdges[i] = edges[i];
	  if (edges[i] == null) return;
	  char reverseColor = MetadataImpl.getReverseEdgeColor(edges[i].getColor());
	  if (reverseColor != edges[i].getColor()){
		  reverseEdges[i] = edges[i].copy();
		  reverseEdges[i].setColor(reverseColor);
	  }
  }

  public void setNumFactors(int numFactors) {
    factors = new int[numFactors];
  }

  public int[] factors() {
    return factors;
  }

  public void addFactors(int[] newFactors) {
    if (newFactors.length != factors.length) {
      throw new RuntimeException("incorrect factor length");
    }
    for (int i=0; i<factors.length; i++) {
      factors[i] += newFactors[i];
    }
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
    edges = new Edge[nodeCount*(nodeCount-1)/2];
    if (useReverseEdges){
    	reverseEdges = new Edge[edges.length];//for multiple site models
    } 
    else {
    	reverseEdges = null;
    }
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
    if (nodes.length == other.nodeCount()) {
      List<Edge> otherEdges = other.edges();
      List<Edge> myEdges = edges();
      if (myEdges.size() == otherEdges.size()) {
        // check store
        int storeOrder = store.compareTo(other.getStore());
        if (storeOrder != 0) {
          return storeOrder;
        }
        // check nodes
        for (byte nodeId = 0; nodeId < nodes.length; nodeId++) {
          int nodeOrder = nodes[nodeId].compareTo(other.getNode(nodeId));
          if (nodeOrder != 0) {
            return nodeOrder;
          }
        }
        // check edges
        for (int iEdge = 0; iEdge < myEdges.size(); iEdge++) {
          Edge edge = myEdges.get(iEdge);
          int edgeOrder = edge.compareTo(otherEdges.get(iEdge));
          if (edgeOrder != 0) {
            return edgeOrder;
          }
        }
        return 0;
      }
      else if (myEdges.size() < otherEdges.size()) {
        return -1;
      }
      else {
        return 1;
      }
    }
    else if (nodes.length < other.nodeCount()) {
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
    Edge[] edgesCopy = new Edge[edges.length];
    for (byte edgeId = 0; edgeId < edges.length; edgeId++) {
      if (edges[edgeId] != null) edgesCopy[edgeId] = edges[edgeId].copy();
    }
    Graph g = new GraphImpl(nodesCopy, edgesCopy, store.copy(), coefficient.copy());
    g.setNumFactors(factors.length);
    g.addFactors(factors);
    return g;
  }

  protected void createEdges() {

    for (byte i = 0; i < nodes.length; i++) {
      for (byte j = (byte) (i + 1); j < nodes.length; j++) {
        byte edgeId = getEdgeId(i, j);
        if (store.testBit(edgeId)) {
          edges[edgeId] = GraphFactory.createEdge(edgeId);
        }
      }
    }
    createReverseEdges();
  }

  public void deleteEdge(byte fromNode, byte toNode) {

    deleteEdge(getEdgeId(fromNode, toNode));
  }

  public void deleteEdge(byte edgeId) {
	edgeId = (byte) (edgeId%edges.length);
    edgeList = null;
    edges[edgeId] = null;
    store.clearBit(edgeId);
  }

  public List<Edge> edges() {
    if (edgeList != null) return edgeList;
    
    edgeList = new ArrayList<Edge>(edgeCount());
    for (byte edgeId=0; edgeId<edges.length; edgeId++) {
      if (edges[edgeId] != null) edgeList.add(edges[edgeId]);
    }
    return edgeList;
  }

  public String edgesToString() {

    // use curly brackets for a set
    String result = "";
    // boolean isMono = (getColors(TYPE_EDGE_ANY).size() == 1);
    boolean first = true;
    for (byte edgeId=0; edgeId<edges.length; edgeId++) {
      if (edges[edgeId] == null) continue;
      if (!first) {
        result += ", ";
      }
      result += "(" + getFromNode(edgeId) + "," + getToNode(edgeId) + ")";
      // if (!isMono) {
      result += edges[edgeId].getColor();
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
      if (edgeCount() != other.edgeCount()) {
        return false;
      }
      // check nodes
      for (int nodeId = 0; nodeId < nodes.length; nodeId++) {
        if (!nodes[nodeId].equals(other.nodes().get(nodeId))) {
          return false;
        }
      }
      // check edges
      for (byte edgeId=0; edgeId<edges.length; edgeId++) {
        if (edges[edgeId] != null && !edges[edgeId].equals(other.getEdge(edgeId))) {
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
    String result = "";
    if (factors.length > 0) {
      result = "/Fac("+factors[0];
      for (int i=1; i<factors.length; i++) {
        result += ","+factors[i];
      }
      result += ")";
    }
    result += "/R";
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
    for (byte edgeId=0; edgeId<edges.length; edgeId++) {
      Edge edge = edges[edgeId];
      if (edge == null) continue;
      char c = edge.getColor();
      if (useReverseEdges) {
        char rc = MetadataImpl.getReverseEdgeColor(c);
        if (c > rc) c = rc;
      }
      Byte count = edgeColorMap.get(c);
      if (count == null) {
        edgeColorMap.put(c, (byte) 1);
      }
      else {
        edgeColorMap.put(c, (byte) (count + 1));
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
    NodeColorVisitor v = new NodeColorVisitor(type);
    visitNodes(v);
    return v.getColors();
  }

  public Edge getEdge(byte fromNode, byte toNode) {
	byte id = getEdgeId(fromNode, toNode);
	if (id < edges.length){
		return edges[id];
	}
    return reverseEdges[id-edges.length];
  }

  public Edge getEdge(byte edgeId) {
	if (edgeId < edges.length){
		return edges[edgeId];
	}
	return reverseEdges[edgeId-edges.length];
  }

  public byte edgeCount() {

    return (byte) store.bitCount();
  }

  // (0,1),(1,2),...,(n-1,0),(0,2),(1,3),...(n-2,0),(n-1,1),...,(0,k),(1,k+1),...,(n-2,k-2),(n-1,k-1)
  // for n odd, k=(n-1)/2
  // for n even, k=(n-2)/2, edge list then also includes (0,n/2),(1,n/2+1),...,(n/2-2,n-2),(n/2-1,n-1)
  public byte getEdgeId(byte fromNode, byte toNode) {
	  boolean reverse = false;//for Wertheim multiple association site

    assert (fromNode != toNode);
    if (fromNode > toNode) {
    	reverse = useReverseEdges;
      byte tmpNode = fromNode;
      fromNode = toNode;
      toNode = tmpNode;
    }
    byte diff = (byte)(toNode - fromNode);
    if (diff > nodes.length/2) {
      // we're jumping from an node near the start to node near the end
      // consider instead going from the node near the end to the node near the beginning
      // (so that diff is smaller)
      diff = (byte)(nodes.length - diff);
      fromNode = toNode;
    }

    // first n edges form the outer ring of edges, etc
    byte id = (byte) ((diff-1)*nodes.length + fromNode);
    if (reverse){
    	id += edges.length;
    }
    return id;
  }

  public byte getFromNode(byte edgeId) {
	boolean reverse = edgeId >= edges.length;
	if (reverse){
		edgeId -= edges.length;
		return getToNode(edgeId);
	}
    byte diff = (byte) (edgeId / nodes.length + 1);
    byte fromNode = (byte) (edgeId - (diff-1)*nodes.length);
    byte toNode = (byte)(fromNode + diff);
    if (toNode < nodes.length) {
      return fromNode;
    }
    // we have something like a from=0, to=(n-1) edge.  our above math gives
    // us from=n-1, to=n, so we need to convert back to from=0
    return (byte)(toNode - nodes.length);
  }

  public Node getNode(byte node) {

    return nodes[node];
  }

  public byte nodeCount() {

    return (byte) nodes.length;
  }

  public byte getOutDegree(byte node) {

    byte result = 0;
    byte nodeCount = (byte)nodes.length;
    for (byte i = 0; i < nodeCount; i++) {
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

  public byte getToNode(byte edgeId) {
	boolean reverse = edgeId >= edges.length;
	if (reverse){
		edgeId -= edges.length;
		return getFromNode(edgeId);
	}
    byte diff = (byte) (edgeId / nodes.length + 1);
    byte fromNode = (byte) (edgeId - (diff-1)*nodes.length);
    byte toNode = (byte)(fromNode + diff);
    if (toNode < nodes.length) {
      return toNode;
    }
    // we have something like a from=0, to=(n-1) edge.  our above math gives
    // us from=n-1, to=n, so we need to take our from node to be our to node.
    return fromNode;
  }

  public boolean hasEdge(byte fromNode, byte toNode) {

    return hasEdge(getEdgeId(fromNode, toNode));
  }

  public boolean hasEdge(byte edgeId) {
    edgeId = (byte) (edgeId%edges.length);
    return store.testBit(edgeId);
  }

  public List<Node> nodes() {
    if (nodeList == null) {
      nodeList = Arrays.asList(nodes);
    }
    return nodeList;
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
    edgeList = null;

    byte edgeId = getEdgeId(fromNode, toNode);
    edgeId = (byte) (edgeId%edges.length);
    store.setBit(edgeId);
    Edge edge = GraphFactory.createEdge(edgeId);
    edges[edgeId] = edge;
    if(useReverseEdges){
    	createReverseEdge(edgeId);
    }
  }

  public void putEdge(byte edgeId) {
    edgeList = null;
    edgeId = (byte) (edgeId%edges.length);
    store.setBit(edgeId);
    edges[edgeId] = GraphFactory.createEdge(edgeId);
    if(useReverseEdges){
    	createReverseEdge(edgeId);
    }
  }

  @Override
  public String toString() {

    String str = coefficient.toString();
    if (factors.length > 0) {
      str += "("+factors[0];
      for (int i=1; i<factors.length; i++) {
        str += ","+factors[i];
      }
      str += ")";
    }
    return str + /* " :: " + getSignature() + */" :: " + nodesToString() + " :: "
        + edgesToString();
  }

  public String toSVG(int dim) {

    int nodeR = dim / 12;
    int oX = dim / 2;
    int oY = dim / 2;
    int graphR = dim / 2 - nodeR;

    // first try 0, then A, then use this fixed one-- should work for 0
    double oA = Math.PI / 2 + Math.PI / nodes.length;
    double angle = 2 * Math.PI / nodes.length;

    double[] x = new double[nodes.length];
    double[] y = new double[nodes.length];

    String svgNodes = "";
    for (int i = 0; i < nodes.length; i++) {
      x[i] = oX + graphR * Math.cos(oA + i * angle);
      y[i] = oY + graphR * Math.sin(oA + i * angle);
      if (nodes.length == 1) {
        x[i] = oX;
        y[i] = oY;
      }
      String svgColor = COLOR_MAP.get(nodes[i].getColor());
      if (svgColor == null) {
        if (COLOR_CODES.indexOf(nodes[i].getColor()) == -1) {
            COLOR_CODES.add(nodes[i].getColor());
        }
        svgColor = COLORS.get(COLOR_CODES.indexOf(nodes[i].getColor()));
        COLOR_MAP.put(nodes[i].getColor(), svgColor);
      }
      if (nodes[i].getType() == TYPE_NODE_ROOT) {
        svgNodes += String.format(
            "<circle style=\"fill: white; stroke: %s\" r=\"%d\" cx=\"%.2f\" cy=\"%.2f\"/>\n",
            svgColor, nodeR, x[i], y[i]);
      }
      else {
        svgNodes += String.format(
            "<circle style=\"fill: %s; stroke: black\" r=\"%d\" cx=\"%.2f\" cy=\"%.2f\"/>\n",
            svgColor, nodeR, x[i], y[i]);
      }
    }
    String svgEdges = "";
    for (Edge e : edges()) {
      byte nodeFrom = getFromNode(e.getId());
      byte nodeTo = getToNode(e.getId());
      String svgColor = COLOR_MAP.get(e.getColor());
      if (svgColor == null) {
        if (COLOR_CODES.indexOf(e.getColor()) == -1) {
            COLOR_CODES.add(e.getColor());
        }
        svgColor = COLORS.get(COLOR_CODES.indexOf(e.getColor()));
        COLOR_MAP.put(e.getColor(), svgColor);
      }
      Integer dashLength = DASH_MAP.get(e.getColor());
      if (dashLength == null) {
        dashLength = 0;
      }
      Integer bondOrderI = BOND_ORDER_MAP.get(e.getColor());
      int bondOrder = bondOrderI == null ? 1 : bondOrderI;
      double odx = y[nodeFrom]-y[nodeTo];
      double ody = x[nodeTo]-x[nodeFrom];
      double norm = Math.sqrt(odx*odx+ody*ody);
      odx /= norm;
      ody /= norm;
      for (int i=0; i<bondOrder; i++) {
        int offset = (2*i-(bondOrder-1))*2;
          
        svgEdges += String.format(
            "<line style=\"stroke: %s;%s\" x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>\n",
            svgColor, (dashLength==0) ? "" : " stroke-dasharray:"+dashLength+","+dashLength, x[nodeFrom]+offset*odx, y[nodeFrom]+offset*ody, x[nodeTo]+offset*odx, y[nodeTo]+offset*ody);
      }
    }
    return svgEdges + svgNodes;
  }

  public void visitEdges(EdgeVisitor visitor) {

    for (byte edgeId=0; edgeId<edges.length; edgeId++) {
      if (edges[edgeId] != null && !visitor.visit(edges[edgeId])) {
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
