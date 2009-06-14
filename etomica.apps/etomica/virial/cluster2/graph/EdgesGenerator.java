package etomica.virial.cluster2.graph;

import java.util.List;

public interface EdgesGenerator extends Tagged {

  public Edges next(List<Edges> edgesList);
}