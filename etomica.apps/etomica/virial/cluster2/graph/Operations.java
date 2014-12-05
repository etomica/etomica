/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.util.List;

import etomica.virial.cluster2.graph.isomorphism.Match;

public class Operations {

  /*
   * PCopy (Positive Copy): GraphSet --> GraphSet 
   * 
   * Given a graph set S, returns a graph set S' with a copy of 
   * every graph in S with the same coefficient.
   */
  public static GraphSet PCopy(final GraphSet gset) {

    assert (gset != null);
    List<Edges> list = GraphFactory.createEdgesList(gset.getNodes().count());
    for (Edges e : gset.getEdgesSet()) {
      list.add(e.copy());
    }
    return GraphFactory.simpleGraphSet(gset.getNodes(), list);
  }

  /*
   * NCopy (Negative Copy): GraphSet --> GraphSet 
   * 
   * Given a graph set S, returns a graph set S' with a copy of 
   * every graph in S with the (additive) inverse coefficient.
   */
  public static GraphSet NCopy(final GraphSet gset) {

    assert (gset != null);
    List<Edges> list = GraphFactory.createEdgesList(gset.getNodes().count());
    for (Edges e : gset.getEdgesSet()) {
      list.add(e.ncopy());
    }
    return GraphFactory.simpleGraphSet(gset.getNodes(), list);
  }

  /*
   * Union (GraphSet Union): GraphSet x GraphSet --> GraphSet 
   * 
   * Given two graph sets S1 and S2, Union(S1, S2) is the bag
   * union of S1 and S2.
   */
  public static GraphSet Union(final GraphSet gset1, final GraphSet gset2) {

    assert (gset1 != null && gset2 != null);
    assert (gset1.getNodes().count() == gset2.getNodes().count());
    GraphSet c1 = PCopy(gset1);
    GraphSet c2 = PCopy(gset2);
    List<Edges> list = GraphFactory.createEdgesList(c1.getNodes().count());
    for (Edges e : c1.getEdgesSet()) {
      list.add(e);
    }
    for (Edges e2 : c2.getEdgesSet()) {
      list.add(e2);
    }
    return GraphFactory.simpleGraphSet(c1.getNodes(), list);
  }

  /*
   * IsoFree (Isomorph-Free GraphSet): GraphSet --> GraphSet 
   * 
   * Given a graph set S, IsoFree(S) returns a graph set S' such that 
   * S' no two graphs g1, g2 in S' are isomorphs of one another. For 
   * every graph g' in S', the coefficient of g' is the sum of the 
   * coefficients of every graph g in S such that g is isomorphic to g'.
   */
  public static GraphSet IsoFree(final GraphSet gset) {

    assert (gset != null);
    GraphSet c1 = PCopy(gset);
    List<Edges> list = GraphFactory.createEdgesList(c1.getNodes().count());
    for (Edges e1 : c1.getEdgesSet()) {
      boolean unique = true;
      for (Edges e2 : list) {
        Graph g1 = GraphFactory.simpleGraph(c1.getNodes(), e1);
        Graph g2 = GraphFactory.simpleGraph(c1.getNodes(), e2);
        if (Match.match(g1, g2)) {
          GraphCoefficient coef1 = e1.getMetadata().getCoefficient();
          GraphCoefficient coef2 = e2.getMetadata().getCoefficient();
          e1.getMetadata().setCoefficient(coef1.add(coef2));
          unique = false;
          break;
        }
      }
      if (unique) {
        list.add(e1);
      }
    }
    return GraphFactory.simpleGraphSet(c1.getNodes(), list);
  }

  /*
   * Delete (GraphSet Deletion): GraphSet x GraphSet --> GraphSet 
   * 
   * Given two graph sets S1 and S2, Delete(S1, S2) is the set
   * resulting from removing every graph in S2 from the set S1.
   */
  public static GraphSet Delete(final GraphSet gset1, final GraphSet gset2) {

    assert (gset1 != null && gset2 != null);
    assert (gset1.getNodes().count() == gset2.getNodes().count());
    List<Edges> list = GraphFactory.createEdgesList(gset1.getNodes().count());
    for (Edges e1 : gset1.getEdgesSet()) {
      boolean unique = true;
      for (Edges e2 : gset2.getEdgesSet()) {
        if (e1.equals(e2)) {
          unique = false;
          break;
        }
      }
      if (unique) {
        list.add(e1);
      }
    }
    return GraphFactory.simpleGraphSet(gset1.getNodes(), list);
  }

  /*
   * Sum (GraphSet Sum): GraphSet x GraphSet --> GraphSet 
   * 
   * Given two graph sets S1 and S2, Sum(S1,S2) = IsoFree(Union(S1,S2)).
   */
  public static GraphSet Sum(final GraphSet gset1, final GraphSet gset2) {

    return IsoFree(Union(gset1, gset2));
  }

  /*
   * Sub (GraphSet Sum): GraphSet x GraphSet --> GraphSet 
   * 
   * Given two graph sets S1 and S2, Sub(S1,S2) = IsoFree(Union(S1,NCopy(S2))).
   */
  public static GraphSet Sub(final GraphSet gset1, final GraphSet gset2) {

    return IsoFree(Union(gset1, NCopy(gset2)));
  }
}