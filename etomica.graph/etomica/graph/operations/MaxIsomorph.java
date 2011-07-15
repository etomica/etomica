package etomica.graph.operations;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.property.Property;

/**
 * Returns a set of graphs with each one corresponding to a graph in the
 * original set.  The graph in the returned set is the "greatest" isomorph of
 * the original graph (according to Graph.compareTo).  The graph in the
 * returned set may just be a copy of the original graph.
 *
 * @author Andrew Schultz
 */
public class MaxIsomorph implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert (params == null || params instanceof MaxIsomorphParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph newGraph = apply(g, (MaxIsomorphParameters)params);
      if (newGraph != null) {
        result.add(newGraph);
      }
    }
    return result;
  }

  public Graph apply(Graph g, MaxIsomorphParameters params) {
    Graph result = params.prop.check(g) ? g.copy() : null;
    if (g.nodeCount() < 2) {
      if (result == null) {
        throw new RuntimeException("no happy graph for "+g);
      }
      return result;
    }
    byte nodeCount = g.nodeCount();
    byte[] labels = new byte[nodeCount];
    for (byte i=0; i<labels.length; i++) {
      labels[i] = i;
    }
    
    Relabel relabel = new Relabel();
    RelabelParameters rp = new RelabelParameters(labels);
    // swaps is the number of times we have swapped node i since the last time we swapped node i-1
    // we'll need to swap nodeCount-1-i times
    // 0 1 2 3 4
    // 0 1 2 4 3
    // 0 1 3 2 4
    // 0 1 3 4 2
    // 0 1 4 2 3
    // 0 1 4 3 2
    // 0 2 1 3 4
    // ...
    // 0 3 1 2 4
    // ...
    // 0 4 1 2 3
    // ...
    // 1 0 2 3 4
    // ...
    // 2 0 1 3 4
    // ...
    // 3 0 1 2 4
    // ...
    // 4 0 1 2 3
    // ...
    byte[] swaps = new byte[nodeCount-1];
    while (true) {
      boolean success = false;
      for (byte iNode = (byte)(nodeCount-2); iNode>-1; iNode--) {
        if (swaps[iNode] < nodeCount-iNode-1) {
          byte tmp = labels[iNode];
          byte minMax = nodeCount;
          byte swapNode = -1;
          // we want to swap with the node at a higher index with the next highest value...
          // we'll have to find it
          for (byte jNode=(byte)(iNode+1); jNode<nodeCount; jNode++) {
            if (labels[jNode] > tmp && labels[jNode] < minMax) {
              minMax = labels[jNode];
              swapNode = jNode;
            }
          }
          if (swapNode > -1) {
            labels[iNode] = labels[swapNode];
            labels[swapNode] = tmp;
            swaps[iNode]++;
            for (byte jNode = (byte)(iNode+1); jNode<(byte)(nodeCount-1); jNode++) {
              swaps[jNode] = 0;
            }
            //sort everything above iNode
            if (iNode+1 < nodeCount-1) {
              Arrays.sort(labels, iNode+1, nodeCount);
            }
            success = true;
            break;
          }
        }
      }
      if (!success) {
        if (result == null) {
          throw new RuntimeException("no happy graph for "+g);
        }
        return result;
      }
      Graph pg = relabel.apply(g, rp);
      if (params.prop.check(pg) && (result == null || pg.compareTo(result) > 0)) {
        result = pg;
      }
    }
  }
  
  /**
   * Parameters class for MaxIsomorph which defines a graph property to select
   * for.  The maximum isomorph with this property will be returned.
   */
  public static class MaxIsomorphParameters implements Parameters {
    public final Property prop;
    public MaxIsomorphParameters(Property prop) {
      this.prop = prop;
    }
  }

  /**
   * Parameters instance which does not filter out any graph property.
   */
  public static final MaxIsomorphParameters PARAM_ALL = new MaxIsomorphParameters(new Property() {
    public boolean check(Graph graph) {return true;}
  });
}
