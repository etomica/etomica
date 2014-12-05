/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;

public abstract class AbstractEdgesFilter implements EdgesFilter {

  private EdgesFilter nextFilter = null;
  private List<String> tags = new ArrayList<String>();

  public AbstractEdgesFilter() {

    computeTags();
  }

  public AbstractEdgesFilter(EdgesFilter filter) {

    chain(filter);
  }

  public boolean accept(Edges edges, List<Edges> edgesList) {

    return doAccept(edges, edgesList) && chainedAccept(edges, edgesList);
  }

  public boolean preAccept(Edges edges, List<Edges> edgesList) {

    return doPreAccept(edges, edgesList) && chainedPreAccept(edges, edgesList);
  }

  public final List<String> getTags() {

    return Collections.unmodifiableList(tags);
  }

  public void chain(EdgesFilter filter) {

    if (nextFilter == null) {
      nextFilter = filter;
    }
    else {
      nextFilter.chain(filter);
    }
    computeTags();
  }

  protected boolean chainedAccept(Edges edges, List<Edges> edgesList) {

    return (getNextFilter() == null || getNextFilter().accept(edges, edgesList));
  }

  protected boolean chainedPreAccept(Edges edges, List<Edges> edgesList) {

    return (getNextFilter() == null || getNextFilter().preAccept(edges, edgesList));
  }

  protected abstract boolean doAccept(Edges edges, List<Edges> edgesList);

  /**
   * By default, doAccept should do all the work. Unless there is an obvious
   * cost saving when using the pre-accept filter before the accept filter,
   * there is really no point in separating the logic in two methods.
   */
  protected boolean doPreAccept(Edges edges, List<Edges> edgesList) {

    return true;
  }

  protected abstract String tag();

  protected final EdgesFilter getNextFilter() {

    return nextFilter;
  }

  protected void computeTags() {

    tags.clear();
    tags.add(tag());
    if (getNextFilter() != null) {
      tags.addAll(getNextFilter().getTags());
    }
  }
}