package etomica.virial.cluster2.graph;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.impl.SimpleEdges;
import etomica.virial.cluster2.graph.impl.SimpleEdgesDecoder;

public class EdgesFactory {

  public static Edges createSimpleEdges(Bitmap edges, EdgesDecoder decoder) {

    return new SimpleEdges(edges, decoder);
  }

  public static Edges createNautyEdges(String edges, EdgesDecoder decoder, double coefficient) {

    return new SimpleEdges(BitmapFactory.getBitmap(edges), decoder, coefficient);
  }

  public static EdgesDecoder createSimpleEdgesDecoder(byte nodeCount) {
    
    return new SimpleEdgesDecoder(nodeCount);
  }
}