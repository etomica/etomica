package etomica.virial.cluster2.mvc;

public interface View {

  /**
   * Initializes the UI, independently of the current state of the model.
   */
  public void initializeUI();

  /**
   * Renders this view visible.
   */
  public void display();
}