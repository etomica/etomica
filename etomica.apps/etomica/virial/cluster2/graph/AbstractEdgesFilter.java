package etomica.virial.cluster2.graph;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public abstract class AbstractEdgesFilter implements EdgesFilter {

  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************

  public AbstractEdgesFilter() {

    computetags();
  }

  public AbstractEdgesFilter(AbstractEdgesFilter filter) {

    setNextFilter(filter);
    computetags();
  }

  // ***********************
  // * PUBLIC METHODS
  // ***********************

  @Override
  public boolean accept(Edges edges) {

    return doAccept(edges) && chainedAccept(edges);
  }

  @Override
  public boolean preAccept(Edges edges) {

    return doPreAccept(edges) && chainedPreAccept(edges);
  }

  @Override
  public final Set<String> getTags() {

    return Collections.unmodifiableSet(tags);
  }

  // ***********************
  // * PROTECTED METHODS
  // ***********************

  protected boolean chainedAccept(Edges edges) {

    return (getNextFilter() == null || getNextFilter().accept(edges));
  }

  protected boolean chainedPreAccept(Edges edges) {

    return (getNextFilter() == null || getNextFilter().preAccept(edges));
  }

  protected abstract boolean doAccept(Edges edges);

  protected abstract boolean doPreAccept(Edges edges);

  protected abstract String tag();

  protected final AbstractEdgesFilter getNextFilter() {

    return nextFilter;
  }

  // ***********************
  // * PRIVATE METHODS
  // ***********************

  private void computetags() {

    tags.add(tag());
    if (getNextFilter() != null) {
      tags.addAll(getNextFilter().getTags());
    }
  }

  private void setNextFilter(AbstractEdgesFilter filter) {

    nextFilter = filter;
  }

  // ***********************
  // * PRIVATE FIELDS
  // ***********************

  private AbstractEdgesFilter nextFilter = null;
  private Set<String> tags = new HashSet<String>();
}