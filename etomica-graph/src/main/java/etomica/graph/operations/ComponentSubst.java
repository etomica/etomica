/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.isomorphism.Match;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.property.NumRootNodes;

/**
 * This class handles the substitution g = A where g is a Graph and A is a set
 * of Graph.  Each graph gArg in Arg is split into components (gArg=c1*c2*c3)
 * and if a component ci "matches" g then we return A*product(cj) where j!=i.
 * 
 * The comparisons are made without regard to root points and when constructing
 * the product to return, we pay attention only to the number of root points in
 * ci (if ci has n root nodes, then we make n root nodes in each graph from A).
 * 
 * @author Andrew Schultz
 */
public class ComponentSubst implements Unary {

  protected final SplitGraph splitGrapher;
  protected final MulFlexible mulFlex;
  protected final MulScalar mulScalar;
  
  public ComponentSubst() {
    splitGrapher = new SplitGraph();
    mulFlex = new MulFlexible();
    mulScalar = new MulScalar();
  }
  
  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      result.addAll(apply(g, (ComponentSubstParameters)params));
    }
    return result;
  }
  
  public Set<Graph> apply(Graph g, ComponentSubstParameters params) {
    GraphList result = new GraphList();
    if (g.nodeCount() < params.gComp.nodeCount()) {
      result.add(g.copy());
      return result;
    }
    // split up the graph into its components
    Set<Graph> split = splitGrapher.apply(g);

    // now recombine them.  if we see an mi1 component, we'll
    // use the substitution above to obtain M graphs
    Set<Graph> gRecombine = new HashSet<Graph>();
    boolean success = false;
    for (Graph gSplit : split) {
      Set<Graph> term;
      Graph gsNoRoot = gSplit.copy();
      for (Node node : gsNoRoot.nodes()) {
        node.setType(Metadata.TYPE_NODE_FIELD);
      }
      if (Match.match(params.gComp, gsNoRoot, true)) {
        // we found an graph component to substitute
        success = true;
        term = mulScalar.apply(params.subst, new MulScalarParameters(gSplit.coefficient()));
        // enforce nRoot_term = nRoot_gPSplit
        // this assumes root node location is unimportant
        int nRoot = NumRootNodes.value(gSplit);
        for (Graph gt : term) {
          for (byte i=0; i<nRoot; i++) {
            gt.getNode(i).setType(Metadata.TYPE_NODE_ROOT);
          }
        }
      }
      else {
        term = new GraphList();
        term.add(gSplit);
      }
      if (gRecombine.isEmpty()) {
        gRecombine.addAll(term);
      }
      else {
        // multiply this term back together with the previous terms
        gRecombine = mulFlex.apply(gRecombine, term, params.mfp);
      }
    }
    if (!success) {
      // we found no mi1 graph, just use the original graph
      result.add(g);
    }
    else {
      result.addAll(gRecombine);
    }
    return result;
  }

  public static class ComponentSubstParameters implements Parameters {
    public final Graph gComp;
    public final Set<Graph> subst;
    public final MulFlexibleParameters mfp;
    public ComponentSubstParameters(Graph gComp, Set<Graph> subst, MulFlexibleParameters mfp) {
      this.gComp = gComp.copy();
      for (Node node : this.gComp.nodes()) {
        if (node.getType() != Metadata.TYPE_NODE_FIELD) {
          throw new RuntimeException("no root nodes here!");
        }
      }
      for (Graph g : subst) {
        if (g.nodeCount() != gComp.nodeCount()) {
          throw new RuntimeException("you're substituting different # of nodes!?!?!");
        }
        for (Node node : g.nodes()) {
          if (node.getType() != Metadata.TYPE_NODE_FIELD) {
            throw new RuntimeException("no root nodes here!");
          }
        }
      }
      this.subst = subst;
      this.mfp = mfp;
    }
  }
}
