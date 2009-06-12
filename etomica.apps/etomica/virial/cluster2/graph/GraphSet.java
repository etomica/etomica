package etomica.virial.cluster2.graph;

import java.util.Set;

/**
 * This interface encodes a family of graphs having a common set of nodes. This
 * is particular useful for representing sets of graphs compactly. For instance,
 * an undirected graph with N nodes may have any of 2^(N(N-1)/2) sets of edges.
 * If we need to represent a large subset of these graphs using one Nodes object
 * per graph, the memory overhead is significant. This interface provides a
 * more conservative approach by sharing the set of nodes across all subgraphs.
 * 
 * @author Demian Lessa
 * 
 */
public interface GraphSet extends Tagged {

  public Set<Edges> getEdgesSet();

  public int getSize();

  public Nodes getNodes();

  public void visitEdgesSet(EdgesSetVisitor visitor);
}