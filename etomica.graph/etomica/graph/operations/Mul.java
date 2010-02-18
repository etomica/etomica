package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import static etomica.graph.model.Metadata.*;

//
///*
//* Mul (Graph Multiplication): Graph x Graph --> Graph
//*/
//public static Graph Mul(final Graph g1, final Graph g2) {
//
//// multiplication: take the union of the nodes
//byte totalNodes = 0;
//byte n1RootCount = 0;
//byte n2RootCount = 0;
//byte totalRootNodes = 0;
//for (byte i = 0; i < n1.count(); i++) {
//  if (n1.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
//    n1RootCount++;
//    totalRootNodes++;
//  }
//  totalNodes++;
//}
//for (byte i = 0; i < n2.count(); i++) {
//  if (n2.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
//    if (i == totalRootNodes) {
//      totalNodes++;
//      totalRootNodes++;
//    }
//    n2RootCount++;
//  }
//  else {
//    totalNodes++;
//  }
//}
//// VERY IMPORTANT: this implementation is still representation dependent, as it
//// explores the fact that root nodes are the first nodes in the representation of a
//// cluster; in order to alleviate this issue, it is necessary to support labels on
//// top of the node representations.
////
//// compute the mapping of field nodes in the new graph to nodes in n1 and n2
//byte fieldIndex = totalRootNodes;
//Map<Byte, Byte> fieldMap = new HashMap<Byte, Byte>();
//for (byte i = 0; i < n1.count(); i++) {
//  if (n1.getAttributes(i).isSameClass(GraphFactory.FIELD_NODE_ATTRIBUTES)) {
//    fieldMap.put(fieldIndex, i);
//    fieldIndex++;
//  }
//}
//for (byte i = 0; i < n2.count(); i++) {
//  if (n2.getAttributes(i).isSameClass(GraphFactory.FIELD_NODE_ATTRIBUTES)) {
//    fieldMap.put(fieldIndex, i);
//    fieldIndex++;
//  }
//}
//char[] rootColors = new char[totalRootNodes];
//char[] fieldColors = new char[totalNodes - totalRootNodes];
//// determine the colors of all root nodes
//for (byte i = 0; i < totalRootNodes; i++) {
//  if (i < n1RootCount) {
//    rootColors[i] = n1.getAttributes(i).getColor();
//  }
//  else {
//    rootColors[i] = n2.getAttributes(i).getColor();
//  }
//}
//// determine the colors of all field nodes
//for (int i = 0; i < totalNodes - totalRootNodes; i++) {
//  if (i < n1.count() - n1RootCount) {
//    fieldColors[i] = n1.getAttributes(fieldMap.get(i)).getColor();
//  }
//  else {
//    fieldColors[i] = n2.getAttributes(fieldMap.get(i)).getColor();
//  }
//}
//// now encode all edges in g1 and g2 to the edges in the new graph
//StringBuffer[] neighbors = new StringBuffer[totalNodes - 1];
//for (int i = 0; i < totalNodes; i++) {
//  neighbors[i] = new StringBuffer(getZeroVector(totalNodes - i - 1));
//}
//for (byte i = 0; i < totalNodes; i++) {
//  for (byte j = (byte) (i + 1); j < totalNodes; j++) {
//    // edges from the first graph
//    if (i < n1.count() && j < n1.count() && g1.getEdges().hasEdge(i, j)) {
//      int fromNode = -1;
//      int toNode = -1;
//      if (n1.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
//        fromNode = i;
//      }
//      else {
//        fromNode = fieldMap.get(i - n1RootCount);
//      }
//      if (n1.getAttributes(j).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
//        toNode = j;
//      }
//      else {
//        toNode = fieldMap.get(j - n1RootCount);
//      }
//      neighbors[fromNode].setCharAt(toNode, '1');
//    }
//    // edges from the second graph
//    else if (i >= n1.count() && j >= n1.count()
//        && g2.getEdges().hasEdge((byte) (i - n1.count()), (byte) (j - n1.count()))) {
//      int fromNode = -1;
//      int toNode = -1;
//      if (n2.getAttributes((byte) (i - n1.count())).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
//        fromNode = i - n1.count();
//      }
//      else {
//        fromNode = fieldMap.get(i - n1.count() - n2RootCount);
//      }
//      if (n2.getAttributes((byte) (j - n1.count())).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
//        toNode = j - n1.count();
//      }
//      else {
//        toNode = fieldMap.get(j - n1.count() - n2RootCount);
//      }
//      neighbors[fromNode].setCharAt(toNode, '1');
//    }
//  }
//}
//// create the the set of nodes
//Nodes nodes = GraphFactory.defaultNodes(fieldColors, rootColors);
//// create the the set of edges
//StringBuffer strRep = new StringBuffer();
//for (int i = 0; i < totalNodes; i++) {
//  strRep.append(neighbors[i]);
//}
//EdgesRepresentationFactory factory = EdgesRepresentationFactory.getFactory(nodes.count());
//Edges edges = GraphFactory.simpleEdges(factory.getRepresentation(BitmapFactory.getBitmap(strRep
//    .toString())));
//Coefficient coef1 = g1.getEdges().getMetadata().getCoefficient();
//Coefficient coef2 = g2.getEdges().getMetadata().getCoefficient();
//Coefficient coef = coef1.multiply(coef2);
//edges.getMetadata().setCoefficient(coef);
//return GraphFactory.simpleGraph(nodes, edges);
//}
public class Mul implements Binary {

  public Graph multiply(Graph left, Graph right) {

    // two root nodes with the same nodeId in left and right must have the same color
    for (byte nodeId = 0; nodeId < left.getNodeCount(); nodeId++) {
      if (nodeId >= right.getNodeCount()) {
        break;
      }
      Node leftNode = left.nodes().get(nodeId);
      Node rightNode = left.nodes().get(nodeId);
      if (leftNode.getType() == TYPE_NODE_ROOT && rightNode.getType() == TYPE_NODE_ROOT) {
        assert (leftNode.isSameColor(rightNode));
      }
    }
    // two edges in left and right must not connect the same labeled root nodes
    // NOTE: expensive check!
//    for (byte i = 0; i < n1.count(); i++) {
//      // root nodes of g1
//      if (n1.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
//        for (byte j = 0; j < n1.count(); j++) {
//          // edges (i,j) between root nodes in g1
//          if (i != j && n1.getAttributes(j).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)
//              && g1.getEdges().hasEdge(i, j)) {
//            // corresponding edges (i,j) between nodes in g2
//            if (i < n2.count() && j < n2.count() && g2.getEdges().hasEdge(i, j)) {
//              // edge (i,j) in g2 must not be between root nodes
//              assert (!n2.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES) || !n2
//                  .getAttributes(j).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES));
//            }
//          }
//        }
//      }
//    }
    return null;
  }

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    Unary pcopy = new PCopy();
    return null;
  }
}