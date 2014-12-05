/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesRepresentationFactory;
import etomica.virial.cluster2.graph.GraphFactory;

public class NaiveEdgesGenerator extends AbstractEdgesGenerator {

  private static final String TAG_COMPLETE = "Naive";
  private Bitmap current;
  private Bitmap maxEdges;
  private EdgesRepresentationFactory factory;

  public NaiveEdgesGenerator(EdgesRepresentationFactory edgesFactory,
      EdgesFilter filter) {

    super(filter);
    factory = edgesFactory;
    if (factory.getNodeCount() == 0) {
      current = BitmapFactory.EMPTY;
      maxEdges = BitmapFactory.EMPTY;
    }
    else if (factory.getNodeCount() == 1) {
      current = BitmapFactory.ZERO;
      maxEdges = BitmapFactory.ZERO;
    }
    else {
      current = BitmapFactory.getBitmap(factory.getMaxEdges(), false);
      maxEdges = BitmapFactory.getBitmap(factory.getMaxEdges(), true);
    }
  }

  @Override
  protected String getTag() {

    return TAG_COMPLETE;
  }

  @Override
  protected Edges push() {

    if (current == BitmapFactory.EMPTY) {
      return null;
    }
    Edges e = GraphFactory.simpleEdges(factory.getRepresentation(current));
    if (current.equals(maxEdges)) {
      current = BitmapFactory.EMPTY;
    }
    else {
      current = current.copy();
      current.inc();
    }
    return e;
  }
}