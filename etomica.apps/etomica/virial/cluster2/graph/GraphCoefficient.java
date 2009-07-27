package etomica.virial.cluster2.graph;

public interface GraphCoefficient {

  public int getValue1();

  public int getValue2();

  public void setValue1(int newValue);

  public void setValue2(int newValue);
  
  public int getSign();

  public GraphCoefficient switchSign();

  public void inc();

  public GraphCoefficient add(GraphCoefficient value);
}
