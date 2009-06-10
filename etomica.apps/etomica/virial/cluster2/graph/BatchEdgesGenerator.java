package etomica.virial.cluster2.graph;

import java.util.Set;

/**
 * This is a simple generator that generates a small number of new edges for
 * every edge passed as an argument to the generate method. An example of a
 * BatchEdgesGenerator could be one that generates the complement of the given
 * edges value.
 * 
 * @author Demian Lessa
 */
public interface BatchEdgesGenerator extends Tagged {

  public Set<Edges> generate(Edges edges);
}