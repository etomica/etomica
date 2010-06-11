package etomica.graph.operations;

public class DeleteEdgeParameters implements Parameters {

  private char color;

  public DeleteEdgeParameters(char color) {

    this.color = color;
  }

  public char color() {

    return color;
  }
}
