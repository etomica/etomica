/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public abstract class CartesianIterator implements GraphIterator {

  private Graph outerGraph = null;
  private GraphIterator outerIterator = null;
  private GraphIterator innerIterator = null;

  protected void bootstrap() {

    setOuterIterator(createOuterIterator());
    if (getOuterIterator().hasNext()) {
      setOuterGraph(getOuterIterator().next());
    }
    setInnerIterator(createInnerIterator());
  }

  // we just computed a graph in the cartesian; it is time to advance to the next
  // positions in the inner and outer iterators
  protected void advance() {

    // if the inner iteration is over, we must advance the outer iteration
    if (!getInnerIterator().hasNext()) {
      setOuterGraph(null);
      // if there is no next outer graph, then we are done;
      // otherwise, we get a new outer graph and a new inner iterator
      if (getOuterIterator().hasNext()) {
        setOuterGraph(getOuterIterator().next());
        setInnerIterator(createInnerIterator());
      }
    }
  }

  protected abstract Graph combineGraphs(Graph outer, Graph inner);

  public abstract GraphIterator createInnerIterator();

  public abstract GraphIterator createOuterIterator();

  public GraphIterator getInnerIterator() {

    return innerIterator;
  }

  public Graph getOuterGraph() {

    return outerGraph;
  }

  public GraphIterator getOuterIterator() {

    return outerIterator;
  }

  // The cartesian product can continue for as long as there exists an outer graph
  // and an inner graph can be obtained from the inner iterator
  public boolean hasNext() {

    return getOuterGraph() != null && getInnerIterator().hasNext();
  }

  public Graph next() {

    if (hasNext()) {
      Graph result = combineGraphs(getOuterGraph(), getInnerIterator().next());
      advance();
      return result;
    }
    return null;
  }

  public void remove() {

    // no-op
  }

  protected void setInnerIterator(GraphIterator iterator) {

    innerIterator = iterator;
  }

  protected void setOuterGraph(Graph graph) {

    outerGraph = graph;
  }

  protected void setOuterIterator(GraphIterator iterator) {

    outerIterator = iterator;
  }
}