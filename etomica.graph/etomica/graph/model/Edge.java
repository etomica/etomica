package etomica.graph.model;

public interface Edge {

  public char getColor();

  public byte getId();

  public Metadata getMetadata();

  public char getType();

  public boolean isCompatible(Edge other);

  public boolean isSameId(Edge other);

  public void setColor(char color);

  public void setType(char type);
}