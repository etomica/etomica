/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.property.IsBiconnected;

/**
 * This operation performs conversion from a set of biconnected graphs with f
 * edges to a set of diagrams of fully connected diagrams with e and f edges.
 * The operation is performed using the procedure as described in
 * 
 * F. H. Ree and W. G. Hoover, J. Chem. Phys. 41, 1635 (1964).
 * 
 * This procedure only works for the full set of biconnected graphs of a given
 * order (or multiple orders), each with only f edges and a coefficient
 * proportional to the number of ways to permute the nodes that yield a
 * different (although equivalent) graph.
 * 
 * @author Andrew Schultz
 */
public class ReeHoover implements Unary {
  
  protected final IsBiconnected isBi = new IsBiconnected();

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert(params instanceof ReeHooverParameters);
    Set<Graph> result = new HashSet<Graph>();
    MulScalar mulScalar = new MulScalar();
    char eBond = ((ReeHooverParameters)params).eBond;
    for (Graph g : argument) {
      int nodeCount = g.nodeCount();
      int maxEdges = nodeCount*(nodeCount-1)/2;
      int coefficient = 1 + calcReeHooverCoefficient(g, 0, (byte)0);
      if (coefficient != 0) {
        MulScalarParameters msp = new MulScalarParameters(coefficient, 1);
        g = mulScalar.apply(g, msp);
        for (byte edgeId = 0; edgeId<maxEdges; edgeId++) {
          if (!g.hasEdge(edgeId)) {
            g.putEdge(edgeId);
            g.getEdge(edgeId).setColor(eBond);
          }
        }
        result.add(g);
      }
    }
    return result;
  }

  public int calcReeHooverCoefficient(Graph g, int numAlreadyRemoved, byte edgeStart) {
    // we're going to assume all edges are f-edges.  If not, you shouldn't be using this class!
    int sum = 0;
    Graph c = g.copy();
    List<Edge> edges = g.edges();
    for (byte iEdge = edgeStart; iEdge < edges.size(); iEdge++) {
      byte edgeId = edges.get(iEdge).getId();
      c.deleteEdge(edgeId);
      if (isBi.check(c)) {
        int numRemoved = numAlreadyRemoved + 1;
        sum += (numRemoved % 2 == 0) ? 1 : -1;
        sum += calcReeHooverCoefficient(c, numRemoved, iEdge);
      }
      c.putEdge(edgeId);
    }
    return sum;
  }
  
  public static class ReeHooverParameters implements Parameters {
    public final char eBond;
    public ReeHooverParameters(char eBond) {
      this.eBond = eBond;
    }
  }
}
