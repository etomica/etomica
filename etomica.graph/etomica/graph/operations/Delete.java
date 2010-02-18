package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

//
//
//  /*
//   * Delete (GraphSet Deletion): GraphSet x GraphSet --> GraphSet
//   *
//   * Given two graph sets S1 and S2, Delete(S1, S2) is the set resulting from removing
//   * every graph in S2 from the set S1.
//   */
//  public static GraphSet Delete(final GraphSet gset1, final GraphSet gset2) {
//
//    assert (gset1 != null && gset2 != null);
//    // TODO: *too weak precondition*
//    assert (gset1.getNodes().count() == gset2.getNodes().count());
//    List<Edges> list = GraphFactory.createEdgesList(gset1.getNodes().count());
//    for (Edges e1 : gset1.getEdgesSet()) {
//      boolean unique = true;
//      for (Edges e2 : gset2.getEdgesSet()) {
//        if (e1.equals(e2)) {
//          unique = false;
//          break;
//        }
//      }
//      if (unique) {
//        list.add(e1);
//      }
//    }
//    return GraphFactory.simpleGraphSet(gset1.getNodes(), list);
//  }
//
public class Delete implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    // TODO: implement equals and hash code for the graphs
    return null;
  }
}