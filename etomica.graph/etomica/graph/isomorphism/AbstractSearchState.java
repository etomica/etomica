package etomica.graph.isomorphism;

import java.util.ArrayList;
import java.util.List;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.Node;

public abstract class AbstractSearchState implements SearchState {

  private Graph firstGraph;
  private Graph secondGraph;
  protected byte core_len;
  protected byte n1, n2;
  protected byte[] core_1, core_2;

  protected AbstractSearchState() {

    // default constructor for use by descendant classes
  }

  public AbstractSearchState(Graph g1, Graph g2) {

    assert (g1 != null);
    assert (g2 != null);
    firstGraph = g1;
    secondGraph = g2;
  }

  public void backTrack() {

    // backtracking does nothing by default
  }

  public byte getCoreLen() {

    return core_len;
  }

  public NodePair[] getCoreSet() {

    ArrayList<NodePair> pairList = new ArrayList<NodePair>();
    for (byte i = 0; i < n1; i++) {
      if (core_1[i] != NULL_NODE) {
        pairList.add(new NodePair(i, core_1[i]));
      }
    }
    return pairList.toArray(new NodePair[] {});
  }

  public List<Edge> getE1() {

    return getG1().edges();
  }

  public Edge getE1(byte fromNodeID, byte toNodeID) {

    return getG1().getEdge(fromNodeID, toNodeID);
  }

  public List<Edge> getE2() {

    return getG2().edges();
  }

  public Edge getE2(byte fromNodeID, byte toNodeID) {

    return getG2().getEdge(fromNodeID, toNodeID);
  }

  public Graph getG1() {

    return firstGraph;
  }

  public Graph getG2() {

    return secondGraph;
  }

  public List<Node> getN1() {

    return getG1().nodes();
  }

  public Node getN1(byte nodeID) {

    return getN1().get(nodeID);
  }

  public List<Node> getN2() {

    return getG2().nodes();
  }

  public Node getN2(byte nodeID) {

    return getN2().get(nodeID);
  }

  protected void setG1(Graph g1) {

    firstGraph = g1;
  }

  protected void setG2(Graph g2) {

    secondGraph = g2;
  }
}