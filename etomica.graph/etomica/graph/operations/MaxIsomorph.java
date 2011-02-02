package etomica.graph.operations;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;

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
    assert (params == null);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph newGraph = apply(g);
      if (newGraph != null) {
        result.add(newGraph);
      }
    }
    return result;
  }

  public Graph apply(Graph g) {
    Graph result = g.copy();
    if (g.nodeCount() < 2) {
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
            if (iNode+1 < nodeCount-1) {
              Arrays.sort(labels, iNode+1, nodeCount);
            }
            success = true;
            break;
          }
        }
      }
      if (!success) {
        return result;
      }
      Graph pg = relabel.apply(g, rp);
      if (pg.compareTo(result) > 0) {
        result = pg;
      }
    }
  }
}
