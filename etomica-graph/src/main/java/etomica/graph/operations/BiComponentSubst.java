/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.isomorphism.Match;
import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.IdenticalGraphFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.property.NumFieldNodes;
import etomica.graph.traversal.BCVisitor;

/**
 * This class handles the substitution g = A where g is a Graph and A is a set
 * of Graph.  Each graph gArg in Arg is split into bicomponents (gArg=c1*c2*c3)
 * and if a component ci "matches" g then we replace the bonding arrangement
 * within ci to match each of the graphs from A.
 * 
 * @author Andrew Schultz
 */
public class BiComponentSubst implements Unary {

  protected final MulFlexible mulFlex;
  protected final MulScalar mulScalar;
  
  public BiComponentSubst() {
    mulFlex = new MulFlexible();
    mulScalar = new MulScalar();
  }
  
  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      result.addAll(apply(g, (BiComponentSubstParameters)params));
    }
    return result;
  }
  
  public Set<Graph> apply(Graph g, BiComponentSubstParameters params) {
    byte n = params.gComp.nodeCount();
    if (g.nodeCount() < n) {
      GraphList result = new GraphList();
      result.add(g.copy());
      return result;
    }
    // split up the graph into its components
    List<List<Byte>> biComps = BCVisitor.getBiComponents(g);
    GraphExtract extractor = new GraphExtract();
    Set<Graph> gRecombine = new HashSet<Graph>();
    gRecombine.add(g);
    List<Byte> array = new ArrayList<Byte>();
    byte[] counters = new byte[n];
    byte[] mappedNodes = new byte[n];
    for (List<Byte> biComp : biComps) {
      if (biComp.size() == n) {
        Graph subGraph = extractor.apply(g, biComp);
        Graph gsNoRoot = subGraph.copy();
        for (Node node : gsNoRoot.nodes()) {
          node.setType(Metadata.TYPE_NODE_FIELD);
        }
        
        if (Match.match(params.gComp, gsNoRoot, true)) {
          // we found a graph component to substitute
          // copy from subst into g
          Set<Graph> newG = new HashSet<Graph>();
          for (Graph gr : gRecombine) {
            for (Graph gs : params.subst) {
              // copy from gs into a copy of gr
              if (params.allPermutations) {
                Set<Graph> permSet = new HashSet<Graph>();
                boolean advanced = true;
                int count = 0;
                while (advanced) {
                  for (byte j=0; j<n; j++) {
                    array.add(biComp.get(j));
                  }
                  advanced = false;
                  
                  for (byte j=0; j<n; j++) {
                    // decorate a point
                    mappedNodes[j] = array.remove(counters[j]);
                    if (!advanced) {
                      counters[j]++;
                      if (counters[j] == n-j) {
                        counters[j] = 0;
                      }
                      else {
                        advanced = true;
                      }
                    }
                  }
                  count++;

                  Graph grNew = mulScalar.apply(gr, new MulScalarParameters(gs.coefficient()));
                  for (byte i=0; i<biComp.size(); i++) {
                    byte ii = mappedNodes[i];
                    for (byte j=0; j<i; j++) {
                      byte jj = mappedNodes[j];
                      if (grNew.hasEdge(jj,ii)) {
                        if (gs.hasEdge(j,i)) {
                          grNew.getEdge(jj,ii).setColor(gs.getEdge(j,i).getColor());
                        }
                        else {
                          grNew.deleteEdge(jj,ii);
                        }
                      }
                      else if (gs.hasEdge(j,i)) {
                        grNew.putEdge(jj,ii);
                        grNew.getEdge(jj,ii).setColor(gs.getEdge(j,i).getColor());
                      }
                    }
                  }
                  permSet.add(grNew);
                }
                newG.addAll(mulScalar.apply(permSet, new MulScalarParameters(1, count)));
              }
              else {
                Graph grNew = mulScalar.apply(gr, new MulScalarParameters(gs.coefficient()));
                for (byte i=0; i<biComp.size(); i++) {
                  byte ii = biComp.get(i);
                  for (byte j=0; j<i; j++) {
                    byte jj = biComp.get(j);
                    if (grNew.hasEdge(jj,ii)) {
                      if (gs.hasEdge(j,i)) {
                        grNew.getEdge(jj,ii).setColor(gs.getEdge(j,i).getColor());
                      }
                      else {
                        grNew.deleteEdge(jj,ii);
                      }
                    }
                    else if (gs.hasEdge(j,i)) {
                      grNew.putEdge(jj,ii);
                      grNew.getEdge(jj,ii).setColor(gs.getEdge(j,i).getColor());
                    }
                  }
                }
                newG.add(grNew);
              }
            }
          }
          gRecombine = newG;
        }
      }
    }
    Set<Graph> result = new HashSet<Graph>();
    IdenticalGraphFilter identicalFilter = new IdenticalGraphFilter(new IteratorWrapper(gRecombine.iterator()));
    Set<Graph> filteredResult = new HashSet<Graph>();
    while (identicalFilter.hasNext()) {
      filteredResult.add(identicalFilter.next());
    }
    if (!params.combineLonelyNodes) {
      return filteredResult;
    }
    for (Graph gr : filteredResult) {
      while (true) {
        byte lonelyNode = -1;
        byte rootNode = -1;
        for (byte i=0; i<gr.nodeCount(); i++) {
          if (lonelyNode == -1 && gr.getOutDegree(i) == 0) {
            lonelyNode = i;
          }
          else if (rootNode == -1 && gr.getNode(i).getType() == Metadata.TYPE_NODE_ROOT) {
            rootNode = i;
          }
        }
        if (lonelyNode == -1 || rootNode == -1) {
          break;
        }
        Graph grNew = GraphFactory.createGraph((byte)(gr.nodeCount()-1));
        for (byte i=0; i<gr.nodeCount(); i++) {
          if (i == lonelyNode) continue;
          byte ii = i;
          if (i > lonelyNode) ii--;
          grNew.getNode(ii).setColor(gr.getNode(i).getColor());
          if (i == rootNode) {
            grNew.getNode(ii).setType(gr.getNode(lonelyNode).getType());
          }
          else {
            grNew.getNode(ii).setType(gr.getNode(i).getType());
          }
          for (byte j=0; j<i; j++) {
            if (!gr.hasEdge(j,i)) continue;
            byte jj = j;
            if (j > lonelyNode) jj--;
            grNew.putEdge(jj,ii);
            grNew.getEdge(jj,ii).setColor(gr.getEdge(j,i).getColor());
          }
        }
        grNew.coefficient().multiply(gr.coefficient());
        gr = grNew;
      }
      result.add(gr);
    }
    return result;
  }
  
  public static class BiComponentSubstParameters implements Parameters {
    public final Graph gComp;
    public final Set<Graph> subst;
    public final boolean allPermutations;
    public final boolean combineLonelyNodes;
    public BiComponentSubstParameters(Graph gComp, Set<Graph> subst, boolean allPermutations, boolean combineLonelyNodes) {
      this.gComp = gComp.copy();
      for (Node node : this.gComp.nodes()) {
        if (node.getType() != Metadata.TYPE_NODE_FIELD) {
          throw new RuntimeException("no root nodes here!");
        }
      }
      int nNodes = gComp.nodeCount();
      for (Graph g : subst) {
        if (NumFieldNodes.value(g) != nNodes) {
          throw new RuntimeException("you're substituting different # of nodes!?!?!");
        }
      }
      this.subst = subst;
      this.allPermutations = allPermutations;
      this.combineLonelyNodes = combineLonelyNodes;
    }
  }

}
