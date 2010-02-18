package etomica.graph.traversal;

public interface TraversalVisitor {

  public boolean visit(byte node, byte status);
}