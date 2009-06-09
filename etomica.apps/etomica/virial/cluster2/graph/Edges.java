package etomica.virial.cluster2.graph;

public interface Edges {

  public static final byte EDGE_COLOR_DEFAULT = 0;
  public static final byte EDGE_COLOR_BLACK = 0;
  public static final byte EDGE_COLOR_RED = 1;

  public int count();

  public int count(byte color);

  public EdgesDecoder getDecoder();

  public EdgesMetadata getMetadata();

  public boolean hasColoredEdge(byte color, int edgeID);

  public boolean hasColoredEdge(byte color, int fromNodeID, int toNodeID);

  public boolean hasEdge(int edgeID);

  public boolean hasEdge(int fromNodeID, int toNodeID);
}