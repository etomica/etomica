package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.AbstractEdgesGenerator;
import etomica.virial.cluster2.graph.EdgesDecoder;
import etomica.virial.cluster2.graph.EdgesFactory;
import etomica.virial.cluster2.graph.EdgesFilter;

public class NaiveEdgesGenerator extends AbstractEdgesGenerator {

  private static final String FLAG_COMPLETE = "Complete";

  private EdgesDecoder decoder;
  private Bitmap current = BitmapFactory.upperTriangleBitmap(getNodeCount(), false);
  private final Bitmap maxEdges = BitmapFactory.upperTriangleBitmap(getNodeCount(), true);

  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************

  public NaiveEdgesGenerator(SimpleGraphSet family, EdgesFilter filter) {

    super(family, filter);
    decoder = new SimpleEdgesDecoder(getNodeCount());
  }

  public NaiveEdgesGenerator(SimpleGraphSet family) {

    this(family, null);
  }

  // ***********************
  // * PROTECTED METHODS
  // ***********************

  @Override
  protected String getTag() {

    return NaiveEdgesGenerator.FLAG_COMPLETE;
  }

  @Override
  protected void pop() {

    if (current == BitmapFactory.EMPTY) {
      setTop(null);
      return;
    }
//    if (isStarted() && getTop() == null) {
//      return;
//    }
    setTop(EdgesFactory.createSimpleEdges(current, decoder));
    if (current.equals(maxEdges)) {
      current = BitmapFactory.EMPTY;
    } else {
      current = current.copy();
      current.inc();
    }
  }
}