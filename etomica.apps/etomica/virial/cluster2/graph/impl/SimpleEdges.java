package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesDecoder;
import etomica.virial.cluster2.graph.EdgesMetadata;

/**
 *  4 bytes instance overhead (experimental: CustomTestCase.testObjectMemoryFootprint()).
 * 12 bytes for data members.
 *  8 bytes of instance overhead for the Bitmap and EdgesDecoder private fields.
 * 12 bytes for data members of BitmapOfLong.
 *  4 bytes for data members of SimpleEdgesDecoder.
 *  0 bytes in EdgesMetadata.
 * ------------------------------
 * 40 bytes per SimpleEdges.
 * 
 * In an array, another 4 bytes are required for a reference to each SimpleEdges 
 * instance. Hence, the asymptotic memory usage for each SimpleEdges instance is
 * about 44 bytes/instance when using BitmapOfLong to encode the graphs.
 * 
 */
public class SimpleEdges implements Edges {

  private Bitmap edges = null;
  private EdgesDecoder edgesDecoder;
  private EdgesMetadata edgesMetadata;

  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************

  public SimpleEdges(Bitmap edgesMap, EdgesDecoder decoder) {

    this(edgesMap, decoder, 1.0);
  }

  public SimpleEdges(Bitmap edgesMap, EdgesDecoder decoder, double coefficient) {

    this(edgesMap, decoder, new SimpleEdgesMetadata(coefficient));
  }

  public SimpleEdges(Bitmap edgesMap, EdgesDecoder decoder, EdgesMetadata metadata) {

    edges = edgesMap;
    edgesDecoder = decoder;
    edgesMetadata = metadata;
  }

  // ***********************
  // * PUBLIC METHODS
  // ***********************

  @Override
  public int count() {

    return getEdges().bitCount();
  }

  @Override
  public int count(byte color) {

    return (color == Edges.EDGE_COLOR_DEFAULT ? count() : 0);
  }

  @Override
  public EdgesDecoder getDecoder() {

    return edgesDecoder;
  }

  @Override
  public EdgesMetadata getMetadata() {

    return edgesMetadata;
  }

  @Override
  public boolean hasColoredEdge(byte color, int edgeID) {

    return (color == Edges.EDGE_COLOR_DEFAULT ? hasEdge(edgeID) : false);
  }

  @Override
  public boolean hasColoredEdge(byte color, int fromNodeID, int toNodeID) {

    return (color == Edges.EDGE_COLOR_DEFAULT ? hasEdge(fromNodeID, toNodeID)
        : false);
  }

  @Override
  public boolean hasEdge(int edgeID) {

    return getEdges().testBit(getDecoder().getBitIndex(edgeID));
  }

  @Override
  public boolean hasEdge(int fromNodeID, int toNodeID) {

    return getEdges().testBit(getDecoder().getBitIndex(fromNodeID, toNodeID));
  }

  @Override
  public int hashCode() {

    return getEdges().hashCode();
  }

  @Override
  public String toString() {

    String result = "";
    Bitmap printSet = getEdges().copy();
    while (printSet.bitCount() > 0) {
      int bit = printSet.hsb();
      result += getDecoder().toString(bit);
      printSet.clearBit(bit);
      if (printSet.bitCount() > 0) {
        result += ", ";
      }
    }
    return "<" + result + ">";
  }

  // ***********************
  // * PROTECTED METHODS
  // ***********************

  protected Bitmap getEdges() {

    return edges;
  }
}