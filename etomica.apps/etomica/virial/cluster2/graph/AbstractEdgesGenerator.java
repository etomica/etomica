package etomica.virial.cluster2.graph;

import java.util.Collections;
import java.util.HashSet;
import java.util.NoSuchElementException;
import java.util.Set;

public abstract class AbstractEdgesGenerator implements EdgesGenerator {

  // ***********************
  // * PRIVATE FIELDS
  // ***********************

  private EdgesFilter edgesFilter = null;
  private Set<String> tags = new HashSet<String>();
  private GraphSet graphSet;
  private boolean started = false;
  private Edges top = null;
  private EdgesGenerator chainedGenerator = null;

  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************

  /**
   * Let G be a chained generator. Then, for every g \in G, we define F(g) to be the 
   * set of graphs generated from g. The generation process is adapted for a chained
   * generator as follows:
   * 
   *   a) top = null, stack = empty stack
   *   b) while G.hasNext()
   *        if stack.empty()
   *          g = G.next()
   *          compute F(g)
   *          stack.push g
   *          for each g' \in F(g) stack.push g'
   *        top = stack.pop
   * 
   * The above implementation is generic even in the absence of a chained generator. 
   * Thus, we generalize the generator to use a stack and push/pop of the values to
   * generate.
   * 
   */
  protected AbstractEdgesGenerator(GraphSet graphs, boolean computeTags, EdgesFilter filter, EdgesGenerator generator) {

    graphSet = graphs;
    edgesFilter = filter;
    chainedGenerator  = generator;
    if (computeTags) {
      computeTags();
    }
  }

  protected AbstractEdgesGenerator(GraphSet graphs, boolean computeTags, EdgesFilter filter) {

    this(graphs, computeTags, filter, null);
  }

  protected AbstractEdgesGenerator(GraphSet graphs, boolean computeTags) {

    this(graphs, computeTags, null, null);
  }

  public AbstractEdgesGenerator(GraphSet graphs, EdgesFilter filter) {

    this(graphs, true, filter, null);
  }

  public AbstractEdgesGenerator(GraphSet graphs) {

    this(graphs, true, null, null);
  }

  // ***********************
  // * PUBLIC METHODS
  // ***********************

  @Override
  public final Set<String> getTags() {

    return Collections.unmodifiableSet(tags);
  }

  @Override
  public final boolean hasNext() {

    if (!isStarted()) {
      filteredPop();
    }
    return (getTop() != null);
  }

  @Override
  public final Edges next() throws NoSuchElementException {

    if (!hasNext()) {
      throw new NoSuchElementException();
    }
    Edges next = getTop();
    filteredPop();
    return next;
  }

  // ***********************
  // * PROTECTED METHODS
  // ***********************

  protected final void filteredPop() {

    pop();
    setStarted(true);
    if (getEdgesFilter() != null) {
      while (getTop() != null && (!getEdgesFilter().preAccept(getTop())
          || !getEdgesFilter().accept(getTop()))) {
        pop();
      }
    }
  }

  protected final EdgesFilter getEdgesFilter() {

    return edgesFilter;
  }

  protected abstract String getTag();

  protected Set<String> getInternalTags() {
    
    return tags;
  }

  protected final EdgesGenerator getChainedGenerator() {

    return chainedGenerator;
  }

  protected final GraphSet getGraphSet() {

    return graphSet;
  }

  protected final byte getNodeCount() {

    return getGraphSet().getNodes().count();
  }

  protected final Edges getTop() {

    return top;
  }

  protected final boolean isStarted() {

    return started;
  }

  protected abstract void pop();

  protected final void setStarted(boolean value) {

    started = value;
  }

  protected final void setTop(Edges edges) {

    top = edges;
  }

  protected void computeTags() {

    if (getChainedGenerator() != null) {
      getInternalTags().addAll(getChainedGenerator().getTags());
    }
    getInternalTags().add(getTag());
    if (getEdgesFilter() != null) {
      getInternalTags().addAll(getEdgesFilter().getTags());
    }
  }
}