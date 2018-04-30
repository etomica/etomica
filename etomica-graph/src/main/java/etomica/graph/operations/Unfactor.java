/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphList;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.property.NumRootNodes;

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
      result.addAll(apply(g, (MulFlexibleParameters)params));
    }
    return result;
  }

  public Set<Graph> apply(Graph g, MulFlexibleParameters params) {
    Set<Graph> splitted = graphSplitter.apply(g);
    Set<Graph> result = null;
    // we might run into problems if the first graphs we see are not
    // superimposable (perhaps no root points) but are both superimposable with
    // a later graph (perhaps having multiple root points).  MulFlex probably
    // should handle that but doesn't.
    // if that does happen, we should superimpose one point from the final multiplication.
    // then, re-split the graph into components and try again to superimpose.  continue as
    // long as at least one superimpose happens
    Set<Graph> unhappy = new HashSet<Graph>();
    int nRoot = 0;
    boolean success = false;
    while (true) {
      success = false;
      for (Graph f : splitted) {
        GraphList set1 = new GraphList();
        set1.add(f);
        if (result == null) {
          result = set1;
          nRoot = NumRootNodes.value(f);
        }
        else {
          int fRoots = NumRootNodes.value(f);
          Set<Graph> newResult = mulFlex.apply(result, set1, params);
          // just need one of the resulting graphs, they shouldn't be meaningfully different
          Graph newg1 = newResult.iterator().next();
          int newRoots = NumRootNodes.value(newg1);
          if (newRoots == nRoot + fRoots) {
            // failed to superimpose
            unhappy.add(f);
          }
          else {
            result = newResult;
            nRoot = newRoots;
            success = true;
          }
        }
      }
      if (unhappy.size() == 0) {
        // we were able to superimpose each pair
        break;
      }
      if (!success) {
        // we we not able to superimpose any pair
        // force multiplication and bail
        for (Graph f : unhappy) {
          GraphList set1 = new GraphList();
          set1.add(f);
          result = mulFlex.apply(result, set1, params);
        }
        break;
      }
      // we failed to superimpose at least 1 component.
      // try again so long as we were successful with at least one pair
      Set<Graph> tmp = splitted;
      splitted = unhappy;
      unhappy = tmp;
      unhappy.clear();
    }
    return result;
  }
}
