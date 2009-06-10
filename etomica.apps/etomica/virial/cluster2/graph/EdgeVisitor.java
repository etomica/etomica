package etomica.virial.cluster2.graph;

public interface EdgeVisitor {

  public boolean visit(int node1, int node2);
}
