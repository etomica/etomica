package etomica.virial.cluster2.mvc;

public interface View {

  /**
   * Present the entire state to the view so it can initialize its internal state.
   *
   */
  public void configure(State state);

  /**
   * Renders this view.
   *
   */
  public void display();
}