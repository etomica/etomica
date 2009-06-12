package etomica.virial.cluster2.graph.impl;

import java.util.Collections;
import java.util.HashSet;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.Stack;

import etomica.virial.cluster2.graph.BatchEdgesGenerator;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesGenerator;

public abstract class AbstractEdgesGenerator implements EdgesGenerator {

  private EdgesFilter edgesFilter = null;
  private Set<String> tags = new HashSet<String>();
  private boolean started = false;
  private BatchEdgesGenerator batchGenerator = null;
  private Stack<Edges> stack = new Stack<Edges>();

  protected AbstractEdgesGenerator(boolean computeTags, EdgesFilter filter,
      BatchEdgesGenerator generator) {

    edgesFilter = filter;
    batchGenerator = generator;
    if (computeTags) {
      computeTags();
    }
  }

  protected AbstractEdgesGenerator(boolean computeTags, EdgesFilter filter) {

    this(computeTags, filter, null);
  }

  protected AbstractEdgesGenerator(boolean computeTags) {

    this(computeTags, null, null);
  }

  public AbstractEdgesGenerator(EdgesFilter filter) {

    this(true, filter, null);
  }

  public AbstractEdgesGenerator() {

    this(true, null, null);
  }

  public final Set<String> getTags() {

    return Collections.unmodifiableSet(tags);
  }

  public final boolean hasNext() {

    if (!isStarted()) {
      smartPush();
    }
    return (!stack.empty());
  }

  public final Edges next() throws NoSuchElementException {

    if (stack.empty()) {
      throw new NoSuchElementException();
    }
    Edges next = stack.pop();
    smartPush();
    return next;
  }

  protected final void smartPush() {

    // do we need to push anything at this time?
    if (!stack.empty()) {
      return;
    }
    // the last pop emptied the stack so push one or more set of edges
    Edges e = push();
    setStarted(true);
    // nothing to push
    if (e == null) {
      return;
    }
    // make sure the set of edges satisfies the chain of filters
    if (getEdgesFilter() != null) {
      while (e != null
          && (!getEdgesFilter().preAccept(e) || !getEdgesFilter().accept(e))) {
        e = push();
      }
      // nothing to push
      if (e == null) {
        return;
      }
    }
    // check whether we need to chain the generation
    if (getBatchGenerator() != null) {
      Set<Edges> candidates = getBatchGenerator().generate(e);
      if (getEdgesFilter() == null) {
        stack.addAll(candidates);
      }
      else {
        for (Edges edgs : candidates) {
          if (getEdgesFilter().preAccept(edgs) && getEdgesFilter().accept(edgs)) {
            stack.push(edgs);
          }
        }
      }
    }
    stack.push(e);
  }

  protected final EdgesFilter getEdgesFilter() {

    return edgesFilter;
  }

  protected abstract String getTag();

  protected Set<String> getInternalTags() {

    return tags;
  }

  protected final BatchEdgesGenerator getBatchGenerator() {

    return batchGenerator;
  }

  protected final boolean isStarted() {

    return started;
  }

  protected abstract Edges push();

  protected final void setStarted(boolean value) {

    started = value;
  }

  protected void computeTags() {

    if (getBatchGenerator() != null) {
      getInternalTags().addAll(getBatchGenerator().getTags());
    }
    getInternalTags().add(getTag());
    if (getEdgesFilter() != null) {
      getInternalTags().addAll(getEdgesFilter().getTags());
    }
  }
}