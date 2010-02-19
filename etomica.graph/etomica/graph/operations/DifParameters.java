package etomica.graph.operations;

public class DifParameters implements Parameters {

  private char color;

  public DifParameters(char color) {

    this.color = color;
  }

  public char color() {

    return color;
  }
}
