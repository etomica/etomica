package etomica.virial.cluster2.graph;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.graph.impl.AbstractBitmapRepresentation;
import etomica.virial.cluster2.graph.impl.AdjacencyMatrixRepresentation;
import etomica.virial.cluster2.graph.impl.UpperTriangleRepresentation;

public abstract class EdgesRepresentationFactory {

  // node count for which the factory will create representations
  private byte nodeCount;

  public EdgesRepresentationFactory(byte nodeCount) {

    this.nodeCount = nodeCount;
  }

  public abstract int getCapacity();

  public byte getNodeCount() {

    return nodeCount;
  }

  public abstract EdgesRepresentation getRepresentation();

  public static EdgesRepresentationFactory getFactory(byte nodeCount) {

    return getFactory(GraphFactory.USE_UPPER_TRIANGLE, nodeCount);
  }

  public static EdgesRepresentationFactory getFactory(boolean useUpperTriangle,
      byte nodeCount) {

    if (useUpperTriangle) {
      return new EdgesRepresentationFactory(nodeCount) {

        @Override
        public int getCapacity() {

          return getNodeCount() * (getNodeCount() - 1) / 2;
        }

        @Override
        public EdgesRepresentation getRepresentation() {

          return new UpperTriangleRepresentation(getNodeCount());
        }
      };
    }
    else {
      return new EdgesRepresentationFactory(nodeCount) {

        @Override
        public int getCapacity() {

          return getNodeCount() * getNodeCount();
        }

        @Override
        public EdgesRepresentation getRepresentation() {

          return new AdjacencyMatrixRepresentation(getNodeCount());
        }
      };
    }
  }

  public EdgesRepresentation getRepresentation(Bitmap current) {

    EdgesRepresentation rep = getRepresentation();
    if (rep instanceof AbstractBitmapRepresentation) {
      ((AbstractBitmapRepresentation) rep).setEdgesBitmap(current);
    }
    return rep;
  }
}