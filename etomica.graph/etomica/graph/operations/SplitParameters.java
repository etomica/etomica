package etomica.graph.operations;

public final class SplitParameters implements Parameters {

  private char edgeColor;
  private char newColor0;
  private char newColor1;

  public SplitParameters(char edgeColor, char newColor0, char newColor1) {

    this.edgeColor = edgeColor;
    this.newColor0 = newColor0;
    this.newColor1 = newColor1;
  }

  public char edgeColor() {

    return edgeColor;
  }

  public char newColor0() {

    return newColor0;
  }

  public char newColor1() {

    return newColor1;
  }
}