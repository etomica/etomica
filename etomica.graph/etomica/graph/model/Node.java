package etomica.graph.model;

public interface Node {

  public Node copy();

  public char getColor();

  public byte getId();

  public Metadata getMetadata();

  public char getType();

  public boolean isCompatible(Node other);

  public boolean isSameColor(Node other);

  public boolean isSameId(Node other);

  public void setColor(char color);

  public void setType(char type);
}