package etomica.virial.cluster2.graph.impl;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;

public abstract class AbstractEdgesFilter implements EdgesFilter {

  private EdgesFilter nextFilter = null;
  private Set<String> tags = new HashSet<String>();

  public AbstractEdgesFilter() {

    computetags();
  }

  public AbstractEdgesFilter(EdgesFilter filter) {

    chain(filter);
  }

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

  @Override
  public void chain(EdgesFilter filter) {

    nextFilter = filter;
    computetags();
  }

  protected boolean chainedAccept(Edges edges) {

    return (getNextFilter() == null || getNextFilter().accept(edges));
  }

  protected boolean chainedPreAccept(Edges edges) {

    return (getNextFilter() == null || getNextFilter().preAccept(edges));
  }

  protected abstract boolean doAccept(Edges edges);

  /**
   * By default, doAccept should do all the work. Unless there is an obvious
   * cost saving when using the pre-accept filter before the accept filter,
   * there is really no point in separating the logic in two methods.
   */
  protected boolean doPreAccept(Edges edges) {

    return true;
  }

  protected abstract String tag();

  protected final EdgesFilter getNextFilter() {

    return nextFilter;
  }

  private void computetags() {

    tags.add(tag());
    if (getNextFilter() != null) {
      tags.addAll(getNextFilter().getTags());
    }
  }
}