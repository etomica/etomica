package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesRepresentation;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.GraphFactory;

public class NaiveEdgesGenerator extends AbstractEdgesGenerator {

  private static final String FLAG_COMPLETE = "Naive";
  private EdgesRepresentation representation;
  private Bitmap current;
  private Bitmap maxEdges;

  public NaiveEdgesGenerator(EdgesRepresentation rep, Bitmap first, Bitmap max,
      EdgesFilter filter) {

    super(filter);
    representation = rep;
    current = first;
    maxEdges = max;
  }

  public NaiveEdgesGenerator(EdgesRepresentation rep, Bitmap first, Bitmap max) {

    this(rep, first, max, null);
  }

  @Override
  protected String getTag() {

    return FLAG_COMPLETE;
  }

  @Override
  protected Edges push() {

    if (current == BitmapFactory.EMPTY) {
      return null;
    }
    Edges e = GraphFactory.simpleEdges(current, representation);
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