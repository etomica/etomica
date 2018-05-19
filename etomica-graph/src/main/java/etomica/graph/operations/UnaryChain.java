/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Coefficient;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphList;

/**
 * Facilitates a chain of unary operations.  This class holds a real operation
 * and a second operation which is performed on the result of the first one.
 * The op given to the constructor is performed and then nextOp.  By making
 * nextOp another UnaryChain, an arbitrary number of operations may be chained
 * together.
 * 
 * @author Andrew Schultz
 */
public class UnaryChain implements Unary {

  public UnaryChain(Unary op, Unary nextOp) {
    this.op = op;
    this.nextOp = nextOp;
    gSet1 = new GraphList();
  }
  
  public Set<Graph> apply(Set<Graph> argument, Parameters chainParams) {
    assert (chainParams instanceof ChainParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      gSet1.add(g);
      result.addAll(op.apply(gSet1, ((ChainParameters)chainParams).params));
      gSet1.clear();
    }
    return result;
  }

  public Graph apply(Graph argument, MulScalarParameters params) {
    Graph result = argument.copy();
    Coefficient c = result.coefficient();
    c.multiply(params.factor());
    return result;
  }
  
  protected final Unary op, nextOp;
  protected final Set<Graph> gSet1;
  
  public static class ChainParameters implements Parameters {
    public final Parameters params;
    public final Parameters nextParams;
    
    public ChainParameters(Parameters params, Parameters nextParams) {
      this.params = params;
      this.nextParams = nextParams;
    }
  }
}
