package etomica.virial.cluster2.graph;

import java.util.List;

public interface EdgesFilter extends Tagged {

  // reserved for simple algorithms- O(1), O(N), etc
  public boolean preAccept(Edges edges, List<Edges> edgesList);

  // reserved for more complex algorithms
  public boolean accept(Edges e, List<Edges> edgesList);

  // chain a filter
  public void chain(EdgesFilter filter);
}