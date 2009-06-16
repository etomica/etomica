package etomica.virial.cluster2.graph;

import java.util.List;

import etomica.virial.cluster2.graph.impl.AbstractEdgesFilter;
import etomica.virial.cluster2.graph.isomorphism.Match;

public class FilterFactory {

  // filter flags
  public static final String TAG_RANGE_FILTER = "EdgesRange";
  public static final String TAG_FALSE_FILTER = "FALSE";
  public static final String TAG_TRUE_FILTER = "TRUE";
  public static final String TAG_CONNECTED = "Connected";
  public static final String TAG_ISOMORPH_FREE = "Isomorph-Free";

  // **************************************************
  //
  // A graph is simply a container for nodes and edges.
  //
  // **************************************************
  public Graph simpleGraph(final Nodes nodes, final Edges edges) {

    return new Graph() {

      public Edges getEdges() {

        return edges;
      }

      public Nodes getNodes() {

        return nodes;
      }
    };
  }

  // **************************************************
  //
  // EdgesFilter instances.
  //
  // **************************************************
  public EdgesFilter trueFilter() {

    return new AbstractEdgesFilter() {

      @Override
      protected boolean doAccept(Edges edges, List<Edges> edgesList) {

        return true;
      }

      @Override
      protected String tag() {

        return TAG_TRUE_FILTER;
      }
    };
  }

  public EdgesFilter falseFilter() {

    return new AbstractEdgesFilter() {

      @Override
      protected boolean doAccept(Edges edges, List<Edges> edgesList) {

        return false;
      }

      @Override
      protected String tag() {

        return TAG_FALSE_FILTER;
      }
    };
  }

  public EdgesFilter connectedFilter(final Nodes nodes) {

    return new AbstractEdgesFilter() {

      @Override
      protected boolean doAccept(Edges edges, List<Edges> edgesList) {

        return Algorithms.isConnected(nodes, edges);
      }

      @Override
      protected String tag() {

        return TAG_CONNECTED;
      }
    };
  }

  public EdgesFilter isomorphismFilter(final Nodes nodes,
      final List<Edges> edgesList) {

    return new AbstractEdgesFilter() {

//      private boolean hardStop = false;
      
      @Override
      protected boolean doAccept(Edges edges, List<Edges> edgesList) {

//        Match.graphs++;
//        if (Match.graphs % 100000 == 0) {
//          System.out.println("===== MARK " + Match.graphs + "=====");
//        }
//        if (hardStop) {
//          return false;
//        }
//        if (edgesList.size() == Match.OPTIMAL_ISMORPHS_COUNT[nodes.count() - 1]) {
//          System.out.println("Optimal upper bound before computing complements.");
//          hardStop = true;
//          return false;
//        }
//        int edgesCount = edges.count();
        Graph g2 = simpleGraph(nodes, edges);
        // if the new graph is isomorphic to some generated graph, reject it
        for (Edges e : edgesList) {
          // save some computation by performing a simpler check before checking
          // for isomorphism
          if (Match.match(simpleGraph(nodes, e), g2)) {
            e.getMetadata()
                .setCoefficient(e.getMetadata().getCoefficient() + 1);
            return false;
          }
        }
//        System.out.println(Match.graphs
//            + "::"
//            + edgesList.size()
//            + "::"
//            + Match.calls
//            + "::"
//            + ((Match.graphs - 1 == Match.oldGraphs) ? 1 : Match.calls
//                / (Match.graphs - Match.oldGraphs)) + " [" + edgesCount + "]");
//        Match.oldGraphs = Match.graphs;
//        Match.calls = 0;
        return true;
      }

      @Override
      protected String tag() {

        return TAG_ISOMORPH_FREE + " (" + Match.DEF_ISOMORPHISM_ALGO + ")";
      }
    };
  }

  public EdgesFilter rangeFilter(int minEdges, int maxEdges) {

    return new RangeFilter(minEdges, maxEdges);
  }

  class RangeFilter extends AbstractEdgesFilter {

    private int maxEdges;
    private int minEdges;

    RangeFilter(int min, int max) {

      minEdges = min;
      maxEdges = max;
    }

    @Override
    protected boolean doAccept(Edges edges, List<Edges> edgesList) {

      int count = edges.count();
      return (count >= minEdges) && (count <= maxEdges);
    }

    @Override
    protected String tag() {

      return TAG_RANGE_FILTER + ":" + minEdges + ":" + maxEdges;
    }
  }
}