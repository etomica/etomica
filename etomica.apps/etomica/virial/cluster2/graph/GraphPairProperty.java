package etomica.virial.cluster2.graph;

/**
 * This interface generalizes algorithms that decide whether graphs G and G'
 * have a particular property P. The property checker algorithm must compare the
 * two graphs to decide whether G and G' have property P. An example of an
 * implementing class is one that checks if G and G' are isomorphic.
 */
public interface GraphPairProperty {

  public boolean check(Nodes nodes, Edges edges1, Edges edges2);
}