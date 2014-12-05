/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;

/*
 * The native graph set file consists of one line for the set of nodes N
 * and one line for each edge set Ei such that Gi(N,Ei) is a graph in the
 * set. The domain for the color values is the natural numbers. The file 
 * is defined as follows:
 * 
 * graph_set_file  = <header> <nodes_record>[1] <edge_set_record>[edge_count]
 * header          = <version> <root_node_count> <field_node_count> <edge_count> new_line
 * nodes_record    = <root_nodes_record> <field_nodes_record> new_line
 * edgeset_record  = <coefficients> <edge_set> new_line
 * 
 * root_nodes_record:
 * 
 * string of length root_node_count where the i-th char in the 
 * list corresponds to the color of the i-th root node in the graph.
 * 
 * field_nodes_record
 * 
 * string of length field_node_count where the i-th char in the 
 * list corresponds to the color of the i-th root node in the graph.
 * 
 * coefficients = <int, int>
 * 
 * a pair of integers representing the two integer values used in the
 * actual computation of the graph's coefficient.
 * 
 * edge_set
 * 
 * upper triangle representation of the edge set modified as follows:
 * every missing edge is represented by the character '0' and every
 * present edge is represented by the character that encodes its color.
 *
 */
public class GraphSetIO {

  private static final String ERROR_INVALID_HEADER_RECORD = "Invalid graph file header";
  private static final String ERROR_INVALID_NODES_RECORD = "Invalid graph file nodes record";
  private static final String ERROR_INVALID_EDGESET_RECORD = "Invalid graph file edge set record";
  private static final String ERROR_INVALID_VERSION = "Invalid graph file version";
  private static String TOKEN_VERSION = "v1.0";
  static String TOKEN_SEPARATOR = ":";

  private static void tokenize(final String line, List<String> tokens) {

    tokens.clear();
    StringTokenizer st = new StringTokenizer(line, TOKEN_SEPARATOR);
    while (st.hasMoreElements()) {
      tokens.add(st.nextToken());
    }
  }

  public static GraphSet readAsText(final BufferedReader reader) {

    assert (reader != null);
    ArrayList<String> tokens = new ArrayList<String>();
    try {
      /*
       * Parse the header
       */
      String header = reader.readLine();
      tokenize(header, tokens);
      if (tokens.size() != 4) {
        throw new RuntimeException(ERROR_INVALID_HEADER_RECORD);
      }
      if (!tokens.get(0).equals(TOKEN_VERSION)) {
        throw new RuntimeException(ERROR_INVALID_VERSION);
      }
      int rootNodeCount = Integer.parseInt(tokens.get(1));
      int fieldNodeCount = Integer.parseInt(tokens.get(2));
      int edgeSetsCount = Integer.parseInt(tokens.get(3));
      /*
       * Parse nodes' colors
       */
      char[] rootColors = new char[rootNodeCount];
      char[] fieldColors = new char[fieldNodeCount];
      String nodeColors = reader.readLine();
      if (nodeColors.length() != rootNodeCount + fieldNodeCount) {
        throw new RuntimeException(ERROR_INVALID_NODES_RECORD);
      }
      for (int i = 0; i < rootNodeCount; i++) {
        rootColors[i] = nodeColors.charAt(i);
      }
      for (int i = 0; i < fieldNodeCount; i++) {
        fieldColors[i] = nodeColors.charAt(rootNodeCount + i);
      }
      /*
       * Parse edge sets
       */
      int edgesCount = (rootNodeCount + fieldNodeCount)
          * (rootNodeCount + fieldNodeCount - 1) / 2;
      int[] edgesCoefficients = new int[2 * edgeSetsCount];
      String[] edgeSets = new String[edgeSetsCount];
      for (int i = 0; i < edgeSetsCount; i++) {
        String edgeSet = reader.readLine();
        tokenize(edgeSet, tokens);
        if (tokens.size() != 3) {
          throw new RuntimeException(ERROR_INVALID_EDGESET_RECORD);
        }
        edgesCoefficients[2 * i] = Integer.parseInt(tokens.get(0));
        edgesCoefficients[2 * i + 1] = Integer.parseInt(tokens.get(1));
        edgeSets[i] = tokens.get(2);
      }
      /*
       * Build the GraphSet
       */
      Nodes nodes = GraphFactory.defaultNodes(fieldColors, rootColors);
      EdgesRepresentationFactory factory = EdgesRepresentationFactory
          .getFactory(true, nodes.count());
      List<Edges> edgesList = GraphFactory.createEdgesList(nodes.count());
      for (int i = 0; i < edgeSetsCount; i++) {
        Bitmap store = BitmapFactory.getBitmap(toBitmap(edgeSets[i]));
        EdgesRepresentation rep = factory.getRepresentation(store);
        edgesList.add(GraphFactory.simpleEdges(rep));
        // TODO: edges color support
      }
      return GraphFactory.simpleGraphSet(nodes, edgesList);
    }
    catch (IOException e) {
      e.printStackTrace();
      throw new RuntimeException(e.getMessage());
    }
  }

  private static String toBitmap(String value) {

    String result = "";
    for (int i = 0; i < value.length(); i++) {
      result += (value.charAt(i) == Edges.EDGE_COLOR_EMPTY ? '0' : '1');
    }
    return result;
  }

  public static void writeAsText(final GraphSet gs, final BufferedWriter writer) {

    assert (gs != null && writer != null);
    try {
      List<Integer> rootNodes = new ArrayList<Integer>();
      List<Integer> fieldNodes = new ArrayList<Integer>();
      NodeAttributes rootAttr = GraphFactory.ROOT_NODE_ATTRIBUTES;
      NodeAttributes fieldAttr = GraphFactory.FIELD_NODE_ATTRIBUTES;
      for (int i = 0; i < gs.getNodes().count(); i++) {
        NodeAttributes na = gs.getNodes().getAttributes(i);
        if (na.isSameClass(rootAttr)) {
          rootNodes.add(i);
        }
        else if (na.isSameClass(fieldAttr)) {
          fieldNodes.add(i);
        }
      }
      /*
       * Write header
       */
      String header = TOKEN_VERSION + TOKEN_SEPARATOR + rootNodes.size()
          + TOKEN_SEPARATOR + fieldNodes.size() + TOKEN_SEPARATOR
          + gs.getEdgesSet().size();
      writer.write(header);
      writer.newLine();
      /*
       * Write nodes record
       */
      String nodesRecord = "";
      for (int i = 0; i < rootNodes.size(); i++) {
        nodesRecord += gs.getNodes().getAttributes(rootNodes.get(i)).getColor();
      }
      for (int i = 0; i < fieldNodes.size(); i++) {
        nodesRecord += gs.getNodes().getAttributes(fieldNodes.get(i))
            .getColor();
      }
      writer.write(nodesRecord);
      writer.newLine();
      /*
       * Write edgeset record
       */
      gs.visitEdgesSet(new WriterVisitor(writer, (byte)(rootNodes.size() + fieldNodes.size())));
      writer.flush();
    }
    catch (IOException e) {
      e.printStackTrace();
      throw new RuntimeException(e.getMessage());
    }
  }
}

class WriterVisitor implements EdgesSetVisitor {

  private byte nodeCount;
  private BufferedWriter writer;

  public WriterVisitor(BufferedWriter bw, byte nodeCount) {

    writer = bw;
    this.nodeCount = nodeCount;
  }

  private String toColoredBitmap(Edges edges) {

    String result = "";
    for (int n1 = 0; n1 < nodeCount; n1++) {
      for (int n2 = n1 + 1; n2 < nodeCount; n2++) {
        result += edges.getAttributes(n1, n2).getColor();
      }
    }
    return result;
  }

  public boolean visit(Edges e) {

    try {
      writer
          .write(String.valueOf(e.getMetadata().getCoefficient().getValue1()));
      writer.write(GraphSetIO.TOKEN_SEPARATOR);
      writer
          .write(String.valueOf(e.getMetadata().getCoefficient().getValue2()));
      writer.write(GraphSetIO.TOKEN_SEPARATOR);
      writer.write(toColoredBitmap(e));
      writer.newLine();
    }
    catch (IOException ioe) {
      ioe.printStackTrace();
      throw new RuntimeException(ioe.getMessage());
    }
    return true;
  }
}