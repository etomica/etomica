package etomica.virial.cluster2.graph.impl;

import java.util.*;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesSetVisitor;
import etomica.virial.cluster2.graph.GraphSet;
import etomica.virial.cluster2.graph.EdgesGenerator;
import etomica.virial.cluster2.graph.GraphException;
import etomica.virial.cluster2.graph.Nodes;

public class SimpleGraphSet implements GraphSet {

  // default allocation space is roughly 7% of maximum
  private static final int DYN_ALLOC_NUM = 1;
  private static final int DYN_ALLOC_DEN = 16;
  private ArrayList<Edges> edges;
  private Nodes nodes;

  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************

  public SimpleGraphSet(byte numNodes) throws GraphException {

    this(numNodes, true);
  }

  public SimpleGraphSet(byte numNodes, boolean allocateFull)
      throws GraphException {

    if (numNodes == 0) {
      throw new GraphException(
          "A GraphSet instance must have at least one node (N=0 not allowed).");
    }
    nodes = new SimpleNodes(numNodes);
    if (allocateFull) {
      edges = new ArrayList<Edges>(edgesInstances(numNodes));
    } else {
      edges = new ArrayList<Edges>(edgesInstancesFraction(numNodes));
    }
  }

  public SimpleGraphSet(byte numNodes, int fracNum, int fracDen)
      throws GraphException {

    if (numNodes == 0) {
      throw new GraphException(
          "A GraphSet instance must have at least one node (N=0 not allowed).");
    }
    nodes = new SimpleNodes(numNodes);
    edges = new ArrayList<Edges>(edgesInstances(numNodes) * fracNum / fracDen);
  }

  protected int edgesInstances(byte numNodes) {

    if (numNodes == 1) {
      return 0x00000001;
    }
    return 0x00000001 << (numNodes * (numNodes - 1) / 2);
  }

  protected int edgesInstancesFraction(byte numNodes) {

    return edgesInstances(numNodes) * SimpleGraphSet.DYN_ALLOC_NUM
        / SimpleGraphSet.DYN_ALLOC_DEN;
  }

  // ***********************
  // * PUBLIC METHODS
  // ***********************

  @Override
  public void generateEdgesSet(EdgesGenerator generator) {

    while (generator.hasNext()) {
      edges.add(generator.next());
    }
  }

  @Override
  public Set<Edges> getEdgesSet() {

    return new HashSet<Edges>(edges);
  }

  @Override
  public int getEdgesSetSize() {

    return edges.size();
  }

  @Override
  public Nodes getNodes() {

    return nodes;
  }

  @Override
  public void replaceEdges(Edges edgesToDrop, Edges edgesToAdd) {

    if (getEdgesSet().remove(edgesToDrop)) {
      getEdgesSet().add(edgesToAdd);
    }
  }

  @Override
  public String toString() {

    String result = "";
    for (int i = 0; i < getEdgesSetSize(); i++) {
      result += i + ": " + edges.get(i).toString() + "\n";
    }
    return result;
  }

  @Override
  public void visitEdgesSet(EdgesSetVisitor visitor) {

    for (int i = 0; (i < getEdgesSetSize()) && visitor.visit(edges.get(i)); i++)
      ;
  }
}