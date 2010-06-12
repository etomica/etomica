package etomica.graph.model;

import java.util.List;

public interface Graph extends Comparable<Graph> {

  public Coefficient coefficient();

  /**
   * factors are considered to be variables that go with the diagram, such as
   * rhoA^2 rhoB^3
   * 
   * The array contains the exponents on each of the factors
   */
  public int[] factors();
  public void addFactors(int[] newFactors);
  public void setNumFactors(int numFactors);

  public Graph copy();

  public void deleteEdge(byte fromNode, byte toNode);

  public byte edgeCount();

  public List<Edge> edges();

  public String edgesToString();

  public Edge getEdge(byte edgeId);

  public Edge getEdge(byte fromNode, byte toNode);

  public byte getEdgeId(byte fromNode, byte toNode);

  public byte getFromNode(byte edge);

  public Node getNode(byte node);

  public byte getOutDegree(byte node);

  public byte getOutNode(byte node, byte index);

  public String getSignature();

  public Bitmap getStore();

  public byte getToNode(byte edge);

  public boolean hasEdge(byte fromNode, byte toNode);

  public byte nodeCount();

  public List<Node> nodes();

  public String nodesToString();

  public void putEdge(byte edgeId);

  public Edge putEdge(byte fromNode, byte toNode);

  public String toSVG(int dim);

  public void visitEdges(EdgeVisitor visitor);

  public void visitNodes(NodeVisitor visitor);
}