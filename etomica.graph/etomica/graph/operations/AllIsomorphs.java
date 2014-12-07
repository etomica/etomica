/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.IdenticalGraphFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.property.Property;

/**
 * Returns all graphs that are isomorphs of the original graph(s).  Optionally,
 * condense graphs that are identical.
 *
 * @author Andrew Schultz
 */
public class AllIsomorphs implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert (params instanceof AllIsomorphsParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      result.addAll(apply(g, (AllIsomorphsParameters)params));
    }
    return result;
  }

  public Set<Graph> apply(Graph g, AllIsomorphsParameters params) {
    Set<Graph> result = new HashSet<Graph>();
    if (params.propertyFilter == null || params.propertyFilter.check(g)) result.add(g.copy());
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
        break;
      }
      Graph newGraph = relabel.apply(g, rp);
      if (params.propertyFilter == null || params.propertyFilter.check(newGraph)) {
        result.add(newGraph);
      }
    }
    int nPermutations = result.size();
    CoefficientImpl multiplier = new CoefficientImpl(1, nPermutations);
    for (Graph gr : result) {
      gr.coefficient().multiply(multiplier);
    }
    if (params.onlyUnique) {
      IdenticalGraphFilter identicalFilter = new IdenticalGraphFilter(new IteratorWrapper(result.iterator()));
      Set<Graph> filteredResult = new HashSet<Graph>();
      while (identicalFilter.hasNext()) {
        filteredResult.add(identicalFilter.next());
      }
      result = filteredResult;
    }
    return result;
  }
  
  public static class AllIsomorphsParameters implements Parameters {
    public boolean onlyUnique;
    public Property propertyFilter;
    public AllIsomorphsParameters(boolean onlyUnique) {
      this(onlyUnique, null);
    }
    public AllIsomorphsParameters(boolean onlyUnique, Property propertyFilter) {
      this.onlyUnique = onlyUnique;
      this.propertyFilter = propertyFilter;
    }
  }
}
