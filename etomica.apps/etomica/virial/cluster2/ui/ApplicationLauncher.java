package etomica.virial.cluster2.ui;

public class ApplicationLauncher {

  /**
   * A simple entry point to run the application. In the future, one might want
   * to use the arguments to the main method to modify the way the UI is
   * constructed.
   */
  public static void main(String args[]) {

    ViewFactory.showView(ViewFactory.VN_VIRIAL);
  }
}