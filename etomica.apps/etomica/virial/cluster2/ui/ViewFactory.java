package etomica.virial.cluster2.ui;

import javax.swing.SwingUtilities;

import etomica.virial.cluster2.ui.wizards.Model;

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
      return new WizardView(VN_WIZARD);
    }
    return null;
  }

  public static View createView(Model current) {

    View view = null;
    // select the view to display based on the model
    return view;
  }
  
  private static class ViewRunnable implements Runnable {

    private String viewName = null;

    public ViewRunnable(final String name) {

      viewName = name;
    }

    public void run() {

      View v = ViewFactory.createView(viewName);
      v.activate();
    }
  }
}