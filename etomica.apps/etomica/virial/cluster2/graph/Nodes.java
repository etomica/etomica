package etomica.virial.cluster2.graph;

public interface Nodes {

  public static final byte NODE_COLOR_DEFAULT = 0;
  public static final byte NODE_COLOR_BLACK = 0;
  public static final byte NODE_COLOR_RED = 1;

  public byte count();

  public byte count(byte color);

  public boolean hasColoredNode(byte color);
}