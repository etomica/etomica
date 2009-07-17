package etomica.virial.cluster2.graph;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/*
 * The native graph set file consists of one line for the set of nodes N
 * and one line for each edge set Ei such that Gi(N,Ei) is a graph in the
 * set. The domain for the color values is the natural numbers. The file 
 * is defined as follows:
 * 
 * graph_set_file = <header> <nodes_record> <edge_set_record>*
 * 
 * nodes_record = <root_nodes_record> <field_nodes_record> new_line
 * 
 * edge_set_record = <coefficient> <edge_set> new_line
 * 
 * root_nodes_record:
 * 
 * the string "R=" followed by a semicolon separated list of integers, 
 * where the i-th integers in the list correspond to the color of the 
 * i-th root node in the graph.
 * 
 * field_nodes_record
 * 
 * the string "F=" followed by a semicolon separated list of integers, 
 * where the i-th integers in the list correspond to the color of the 
 * i-th field node in the graph.
 * 
 * coefficient
 * 
 * the coefficient for the graph Gi(N,Ei), given as a pair of integers
 * <int1, int2> where the actual coefficient is int1/int2.
 * 
 * edge_set
 * 
 * a semicolon separated list of integers, where the i-th integers in 
 * the list correspond to the color of the i-th edge in the graph. We
 * only include entries for the upper triangle of the adjacency matrix
 * of the corresponding graph. So, if N is the number of nodes, the
 * total number of values in the list is (N-1)(N-2)/2. The value ZERO
 * indicates the absence of an edge and, since ZERO is not part of the
 * color domain, no ambiguity can arise from this convention.
 *
 */
public interface GraphSetIO {

  public GraphSet readAsNative(final BufferedReader reader);

  public void writeAsNative(final GraphSet gs, final BufferedWriter writer);

  public void writeAsText(final GraphSet gs, final BufferedWriter writer);
}