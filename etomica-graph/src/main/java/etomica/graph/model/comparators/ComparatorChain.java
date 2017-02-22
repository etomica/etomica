/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.comparators;

import java.util.ArrayList;
import java.util.Comparator;

import etomica.graph.model.Graph;

/**
 * Comparator class that can chain together a list of comparators.  The first
 * comparator that recognizes one graph as coming before another is used.  If
 * none of the comparators gives precedence to one of the graphs, this class
 * falls back to the Graph's internal comparison, and then the hashCode.
 * 
 * @author Andrew Schultz
 */
public class ComparatorChain implements Comparator<Graph> {

  public ComparatorChain() {
    comparators = new ArrayList<Comparator<Graph>>();
    comparatorDirections = new ArrayList<Integer>();
  }

  /**
   * Add a comparator that will be used to sort sort graphs.  Comparators
   * will be used in the order in which they are added; later comparators
   * will be used only if earlier comparators return 0.
   */
  public void addComparator(Comparator<Graph> comparator) {
    comparators.add(comparator);
    comparatorDirections.add(1);
  }

  /**
   * Add a comparator that will be used to sort sort graphs in reverse order.
   */
  public void addReverseComparator(Comparator<Graph> comparator) {
    comparators.add(comparator);
    comparatorDirections.add(-1);
  }

  public int compare(Graph g1, Graph g2) {
    for (int i=0; i<comparators.size(); i++) {
      int c = comparators.get(i).compare(g1, g2) * comparatorDirections.get(i);
      if (c != 0) {
        return c;
      }
    }
    // fall back on internal comparison
    // negate what we get here -- isomorphism prefers high-score graphs.
    // we want those same graphs to come first
    int c = -g1.compareTo(g2);
    if (c != 0) return c;
    // graphs are indistinguishable, just return something consistent
    return g1.hashCode() - g2.hashCode();
  }

  protected final ArrayList<Comparator<Graph>> comparators;
  protected final ArrayList<Integer> comparatorDirections;
}