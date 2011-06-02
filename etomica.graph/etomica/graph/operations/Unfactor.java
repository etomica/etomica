package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;

/**
 * This operation will do the reverse of the Factor operation.  Given a graph
 * with root points that is disconnected, it will return a graph that has these
 * root points superimposed with other points, such that the graph is connected
 * if possible.  This is done by first splitting the graph into components and
 * then multiplying them together.
 */
public class Unfactor implements Unary {

  protected final SplitGraph graphSplitter = new SplitGraph();
  protected final MulFlexible mulFlex = new MulFlexible();
  
  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert(params instanceof MulFlexibleParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      result.add(apply(g, (MulFlexibleParameters)params));
    }
    return result;
  }

  public Graph apply(Graph g, MulFlexibleParameters params) {
    Set<Graph> splitted = graphSplitter.apply(g);
    Graph result = null;
    // we might run into problems if the first graphs we see are
    // superimposable (perhaps no root points) but are both superimposable with
    // a later graph (perhaps having multiple root points).  MulFlex probably
    // should handle that but doesn't.
    for (Graph f : splitted) {
      if (result == null) {
        result = f;
      }
      else {
        result = mulFlex.apply(result, f, params);
      }
    }
    return result;
  }
}
