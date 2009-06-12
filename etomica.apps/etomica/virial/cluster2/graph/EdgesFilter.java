package etomica.virial.cluster2.graph;

public interface EdgesFilter extends Tagged {

  // reserved for more complex algorithms
  public boolean accept(Edges edges);

  // reserved for simple algorithms- O(1), O(N), etc
  public boolean preAccept(Edges edges);

  // chain a filter
  public void chain(EdgesFilter filter);
}