package etomica.graph.engine;

public interface Console {

  public void clear();

  public void setReader(ConsoleReader reader);

  public void writeLn();

  public void write(Exception e);

  public void write(String value);

  public void updateBegin();

  public void updateDone();
}