package etomica.virial.cluster2.mvc.launchers;

import javax.swing.SwingUtilities;

import etomica.virial.cluster2.mvc.View;
import etomica.virial.cluster2.mvc.view.ApplicationUI;
import etomica.virial.cluster2.mvc.view.ApplicationView;
import etomica.virial.cluster2.mvc.view.ClusterWizard;

public class ViewFactory {

  public static String VN_VIRIAL = "virialApplication";
  public static String VN_WIZARD = "wizardTest";

  public static void showView(final String viewName) {

    ApplicationUI.configure();
    SwingUtilities.invokeLater(new ViewRunnable(viewName));
  }

  public static View createView(String viewName) {

    if (viewName.equalsIgnoreCase(VN_VIRIAL)) {
      return new ApplicationView();
    }
    if (viewName.equalsIgnoreCase(VN_WIZARD)) {
      return new ClusterWizard(VN_WIZARD);
    }
    return null;
  }

  private static class ViewRunnable implements Runnable {

    private String viewName = null;

    public ViewRunnable(final String name) {

      viewName = name;
    }

    public void run() {

      View v = ViewFactory.createView(viewName);
      // TODO: provide a real state
      v.configure(null);
      v.display();
    }
  }
}