/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.util.ArrayList;
import java.util.List;

import etomica.virial.cluster2.graph.impl.NaiveEdgesGenerator;
import etomica.virial.cluster2.graph.impl.NautyEdgesGenerator;
import etomica.virial.cluster2.graph.impl.SimpleEdgeAttributes;
import etomica.virial.cluster2.graph.impl.SimpleEdges;
import etomica.virial.cluster2.graph.impl.SimpleEdgesMetadata;
import etomica.virial.cluster2.graph.impl.SimpleGraphCoefficient;
import etomica.virial.cluster2.graph.impl.SimpleGraphSet;
import etomica.virial.cluster2.graph.impl.SimpleNodes;
import etomica.virial.cluster2.graph.impl.StoredEdgesGenerator;
import etomica.virial.cluster2.nauty.NautyInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;

public class GraphFactory {

  // default representation is adjacency matrix
  public static boolean USE_UPPER_TRIANGLE = true;
  public static final NodeAttributes FIELD_NODE_ATTRIBUTES = SimpleNodes.fieldNodeAttributes(Nodes.NODE_COLOR_DEFAULT);
  public static final NodeAttributes ROOT_NODE_ATTRIBUTES = SimpleNodes.rootNodeAttributes(Nodes.NODE_COLOR_DEFAULT);
  // a class of edge attributes that are always compatible
  public static final EdgeAttributes DEFAULT_EDGE_ATTRIBUTES = SimpleEdgeAttributes.getAttributes(Edges.EDGE_COLOR_DEFAULT);
  public static final EdgeAttributes EMPTY_EDGE_ATTRIBUTES = SimpleEdgeAttributes.getAttributes(Edges.EDGE_COLOR_EMPTY);
  // default allocation space is roughly 3% of maximum (1/32)
  public static int DYN_ALLOC_NUM = 1;
  public static int DYN_ALLOC_DEN = 32;
  public static boolean DYNAMIC_ALLOCATION = true;
  // default coefficient for a graph
  public static GraphCoefficient DEFAULT_GRAPH_COEFFICIENT = new SimpleGraphCoefficient(1);

  // **************************************************
  //
  // A graph is simply a container for nodes and edges,
  // and a GraphSet is a container for nodes and a set
  // of edges.
  //
  // **************************************************
  public static Graph simpleGraph(final Nodes nodes, final Edges edges) {

    return new Graph() {

      public Edges getEdges() {

        return edges;
      }

      public Nodes getNodes() {

        return nodes;
      }
    };
  }

  public static GraphSet simpleGraphSet(final Nodes nodes,
      final List<Edges> list) {

    return new SimpleGraphSet(nodes, list);
  }

  // **************************************************
  //
  // Edges related methods.
  //
  // **************************************************
  private static EdgesMetadata simpleEdgesMetadata(GraphCoefficient coefficient) {

    return new SimpleEdgesMetadata(coefficient);
  }

  public static Edges simpleEdges(final EdgesRepresentation rep) {

    return new SimpleEdges(rep, simpleEdgesMetadata(DEFAULT_GRAPH_COEFFICIENT));
  }

  /**
   * Each graph string output by the modified nauty implementation is either an
   * upper triangle or an adjacency matrix representation of the graph. Hence,
   * the client code must make sure the flags passed to nauty are consistent
   * with the flags used to create the representation object.
   */
  public static Edges nautyEdges(final EdgesRepresentation rep,
      GraphCoefficient coefficient) {

    return new SimpleEdges(rep, simpleEdgesMetadata(coefficient));
  }

  public static Edges complementEdges(final Edges edges) {

    return new SimpleEdges(edges.getRepresentation().complement(),
        simpleEdgesMetadata(edges.getMetadata().getCoefficient()));
  };

  // **************************************************
  //
  // Nodes related methods.
  //
  // **************************************************
  public static Nodes defaultNodes(byte nodeCount) {

    return new SimpleNodes(nodeCount, (byte) 0);
  }

  public static Nodes defaultNodes(byte fieldNodes, byte rootNodes) {

    return new SimpleNodes(fieldNodes, rootNodes);
  }

  public static Nodes defaultNodes(char[] fieldColors, char[] rootColors) {

    return new SimpleNodes(fieldColors, rootColors);
  }

  public static Nodes emptyNodes() {

    return new SimpleNodes((byte) 0, (byte) 0);
  }

  // **************************************************
  //
  // GraphFamily helper methods.
  //
  // **************************************************
  private static int maxInstances(byte numNodes) {

    if (numNodes == 1) {
      return 0x00000001;
    }
    return 0x00000001 << (numNodes * (numNodes - 1) / 2);
  }

  public static List<Edges> createEdgesList(byte numNodes) {

    int listSize = maxInstances(numNodes);
    if (DYNAMIC_ALLOCATION) {
      listSize = listSize * DYN_ALLOC_NUM / DYN_ALLOC_DEN;
    }
    return new ArrayList<Edges>(listSize);
  }

  // **************************************************
  //
  // GraphFamily instances generated naively.
  //
  // **************************************************
  public static GraphSet naiveGraphSet(final Nodes nodes,
      final EdgesFilter filter) {

    return naiveGraphSet(nodes, filter, false);
  }

  public static GraphSet naiveGraphSet(final Nodes nodes,
      final EdgesFilter filter, boolean isomorphFree) {

    EdgesRepresentationFactory factory = EdgesRepresentationFactory
        .getFactory(nodes.count());
    List<Edges> edgesList = createEdgesList(nodes.count());
    EdgesFilter ef = filter;
    if (isomorphFree) {
      EdgesFilter isof = (new FilterFactory()).isomorphismFilter(nodes,
          edgesList);
      if (ef != null) {
        ef.chain(isof);
      }
      else {
        ef = isof;
      }
    }
    EdgesGenerator generator = new NaiveEdgesGenerator(factory, ef);
    int i = 0;
    Edges e = generator.next(edgesList);
    while (e != null) {
      i++;
      edgesList.add(e);
      e = generator.next(edgesList);
    }
    SimpleGraphSet result = (SimpleGraphSet)simpleGraphSet(nodes, edgesList);
    result.setTags(generator.getTags());
    return result;
  }

  // **************************************************
  //
  // GraphFamily instances generated from stored nauty
  // graphs in graph6 format.
  //
  // **************************************************
  public static GraphSet storedGraphSet(final Nodes nodes,
      final EdgesFilter filter) {

    return storedGraphSet(nodes, filter, false);
  }

  public static GraphSet storedGraphSet(final Nodes nodes,
      final EdgesFilter filter, boolean isomorphFree) {

    EdgesRepresentationFactory factory = EdgesRepresentationFactory
        .getFactory(nodes.count());
    List<Edges> edgesList = createEdgesList(nodes.count());
    EdgesFilter ef = filter;
    if (isomorphFree) {
      EdgesFilter isof = (new FilterFactory()).isomorphismFilter(nodes,
          edgesList);
      if (ef != null) {
        ef.chain(isof);
      }
      else {
        ef = isof;
      }
    }
    EdgesGenerator generator = new StoredEdgesGenerator(factory, ef);
    int i = 0;
    Edges e = generator.next(edgesList);
    while (e != null) {
      i++;
      edgesList.add(e);
      e = generator.next(edgesList);
    }
    SimpleGraphSet result = (SimpleGraphSet)simpleGraphSet(nodes, edgesList);
    result.setTags(generator.getTags());
    return result;
  }

  // **************************************************
  //
  // GraphFamily instances generated using nauty.
  //
  // **************************************************
  public static GraphSet nautyGraphSet(final Nodes nodes,
      final EdgesFilter filter, final ProcessWrapper nauty) {

    EdgesRepresentationFactory factory = EdgesRepresentationFactory.getFactory(
        ((NautyInfo) nauty.getProcessInfo()).isUpperTriangle(), nodes.count());
    EdgesGenerator generator = new NautyEdgesGenerator(factory, nauty, filter);
    List<Edges> edgesList = createEdgesList(nodes.count());
    Edges e = generator.next(edgesList);
    while (e != null) {
      edgesList.add(e);
      e = generator.next(edgesList);
    }
    SimpleGraphSet result = (SimpleGraphSet)simpleGraphSet(nodes, edgesList);
    result.setTags(generator.getTags());
    return result;
  }

  public static GraphCoefficient defaultCoefficient(int coefficient) {

    return new SimpleGraphCoefficient(coefficient);
  }
}