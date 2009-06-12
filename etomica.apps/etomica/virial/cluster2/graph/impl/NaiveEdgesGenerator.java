package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesRepresentationFactory;
import etomica.virial.cluster2.graph.GraphFactory;

public class NaiveEdgesGenerator extends AbstractEdgesGenerator {

  private static final String FLAG_COMPLETE = "Naive";
  private EdgesRepresentationFactory factory;
  private Bitmap current;
  private Bitmap maxEdges;

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
      current = BitmapFactory.getBitmap(factory.getCapacity(), false);
      maxEdges = BitmapFactory.getBitmap(factory.getCapacity(), true);
    }
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