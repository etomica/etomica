/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.util.List;

import etomica.virial.cluster2.graph.impl.AbstractEdgesFilter;
import etomica.virial.cluster2.graph.isomorphism.Match;

public class FilterFactory {

	// filter flags
	public static final String TAG_RANGE_FILTER = "Edges Range";
	public static final String TAG_FALSE_FILTER = "FALSE";
	public static final String TAG_TRUE_FILTER = "TRUE";
	public static final String TAG_BICONNECTED = "Biconnected";
	public static final String TAG_NODAL_POINT = "Nodal Point";
  public static final String TAG_ARTICULATION_POINT = "Articulation Point";
	public static final String TAG_CONNECTED = "Connected";
	public static final String TAG_ISOMORPH_FREE = "Isomorph-Free";
	public static final String TAG_NO_ROOT_EDGES = "No Root Edges";
	public static final String TAG_ISOMORPH_FIELD_RANGE_FILTER = "Field Edges Range";
	public static final String TAG_ISOMORPH_ROOT_RANGE_FILTER = "Root Edges Range";

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

	public EdgesFilter biconnectedFilter(final Nodes nodes) {

		return new AbstractEdgesFilter() {

			@Override
			protected boolean doAccept(Edges edges, List<Edges> edgesList) {

				return Algorithms.isBiconnected(nodes, edges);
			}

			@Override
			protected String tag() {

				return TAG_BICONNECTED;
			}
		};
	}

	public EdgesFilter nodalPointFilter(final Nodes nodes) {

		return new AbstractEdgesFilter() {

			@Override
			protected boolean doAccept(Edges edges, List<Edges> edgesList) {

				return !Algorithms.hasNodalPoint(nodes, edges);
			}

			@Override
			protected String tag() {

				return TAG_NODAL_POINT;
			}
		};
	}

  public EdgesFilter articulationPointFilter(final Nodes nodes) {

    return new AbstractEdgesFilter() {

      @Override
      protected boolean doAccept(Edges edges, List<Edges> edgesList) {

        return !Algorithms.hasArticulationPoint(nodes, edges);
      }

      @Override
      protected String tag() {

        return TAG_ARTICULATION_POINT;
      }
    };
  }

  public EdgesFilter isomorphismFilter(final Nodes nodes,
			final List<Edges> edgesList) {

		return new AbstractEdgesFilter() {

			@Override
			protected boolean doAccept(Edges edges, List<Edges> edgesList) {

				Match.graphs++;
				if (Match.graphs % 100000 == 0) {
					System.out.println("===== MARK " + Match.graphs + "=====");
				}
				// int edgesCount = edges.count();
				Graph g2 = GraphFactory.simpleGraph(nodes, edges);
				// if the new graph is isomorphic to some generated graph,
				// reject it
				for (Edges e : edgesList) {
					// save some computation by performing a simpler check
					// before checking
					// for isomorphism
					if (Match.match(GraphFactory.simpleGraph(nodes, e), g2)) {
						e.getMetadata().getCoefficient().inc();
						return false;
					}
				}
				// System.out.println(Match.graphs
				// + "::"
				// + edgesList.size()
				// + "::"
				// + Match.calls
				// + "::"
				// + ((Match.graphs - 1 == Match.oldGraphs) ? 1 : Match.calls
				// / (Match.graphs - Match.oldGraphs)) + " [" + edgesCount +
				// "]");
				// Match.oldGraphs = Match.graphs;
				// Match.calls = 0;
				return true;
			}

			@Override
			protected String tag() {

				return TAG_ISOMORPH_FREE + " (" + Match.DEF_ISOMORPHISM_ALGO
						+ ")";
			}
		};
	}

	public EdgesFilter rangeFilter(int minEdges, int maxEdges) {

		return new RangeFilter(minEdges, maxEdges);
	}

	public EdgesFilter rangedFieldEdgesFilter(Nodes nodes, int minEdges,
			int maxEdges) {

		return new RangedFieldEdgesFilter(nodes, minEdges, maxEdges);
	}

	public EdgesFilter rangedRootEdgesFilter(Nodes nodes, int minEdges,
			int maxEdges) {

		return new RangedRootEdgesFilter(nodes, minEdges, maxEdges);
	}

	public EdgesFilter rootEdgesFilter(final Nodes nodes) {

		return new AbstractEdgesFilter() {

			@Override
			protected boolean doAccept(Edges edges, List<Edges> edgesList) {

				for (int i = 0; i < nodes.count(); i++) {
					for (int j = i; j < nodes.count(); j++) {
						if ((i != j)
								&& nodes.getAttributes(i).isSameClass(
										GraphFactory.ROOT_NODE_ATTRIBUTES)
								&& nodes.getAttributes(j).isSameClass(
										GraphFactory.ROOT_NODE_ATTRIBUTES)
								&& edges.hasEdge(i, j)) {
							return false;
						}
					}
				}
				return true;
			}

			@Override
			protected String tag() {

				return TAG_NO_ROOT_EDGES;
			}
		};
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

	class RangedFieldEdgesFilter extends AbstractEdgesFilter {

		private int maxEdges;
		private int minEdges;
		private Nodes nodes;

		RangedFieldEdgesFilter(Nodes nodes, int min, int max) {

			minEdges = min;
			maxEdges = max;
			this.nodes = nodes;
		}

		@Override
		protected boolean doAccept(Edges edges, List<Edges> edgesList) {

			int count = 0;
			for (int i = 0; i < nodes.count(); i++) {
				for (int j = i; j < nodes.count(); j++) {
					if ((i != j)
							&& nodes.getAttributes(i).isSameClass(
									GraphFactory.FIELD_NODE_ATTRIBUTES)
							&& nodes.getAttributes(j).isSameClass(
									GraphFactory.FIELD_NODE_ATTRIBUTES)
							&& edges.hasEdge(i, j)) {
						count++;
					}
				}
			}
			return (count >= minEdges) && (count <= maxEdges);
		}

		@Override
		protected String tag() {

			return TAG_ISOMORPH_FIELD_RANGE_FILTER + ":" + minEdges + ":"
					+ maxEdges;
		}
	}

	class RangedRootEdgesFilter extends AbstractEdgesFilter {

		private int maxEdges;
		private int minEdges;
		private Nodes nodes;

		RangedRootEdgesFilter(Nodes nodes, int min, int max) {

			minEdges = min;
			maxEdges = max;
			this.nodes = nodes;
		}

		@Override
		protected boolean doAccept(Edges edges, List<Edges> edgesList) {

			int count = 0;
			for (int i = 0; i < nodes.count(); i++) {
				for (int j = i; j < nodes.count(); j++) {
					if ((i != j)
							&& (nodes.getAttributes(i).isSameClass(
									GraphFactory.ROOT_NODE_ATTRIBUTES) || nodes
									.getAttributes(j).isSameClass(
											GraphFactory.ROOT_NODE_ATTRIBUTES))
							&& edges.hasEdge(i, j)) {
						count++;
					}
				}
			}
			return (count >= minEdges) && (count <= maxEdges);
		}

		@Override
		protected String tag() {

			return TAG_ISOMORPH_ROOT_RANGE_FILTER + ":" + minEdges + ":"
					+ maxEdges;
		}
	}
}