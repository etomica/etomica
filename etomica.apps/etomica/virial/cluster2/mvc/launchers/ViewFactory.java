/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.launchers;

import javax.swing.SwingUtilities;

import etomica.virial.cluster2.mvc.View;
import etomica.virial.cluster2.mvc.view.ApplicationUI;
import etomica.virial.cluster2.mvc.view.ApplicationView;
import etomica.virial.cluster2.mvc.view.ClusterWizard;
import etomica.virial.cluster2.mvc.view.ClusterWizardController;

public class ViewFactory {

  public static String VN_VIRIAL = "virialApplication";
  public static String VN_EMPTY_WIZARD = "Wizard Demo";
  public static String VN_CLUSTER_WIZARD = "Cluster Wizard";

  public static void showView(final String viewName) {

    ApplicationUI.configure();
    // the wizard controller implements its own runnable, so no need to wrap here
    if (viewName.equalsIgnoreCase(VN_CLUSTER_WIZARD)) {
      SwingUtilities.invokeLater(new ClusterWizardController());
    }
    else {
      SwingUtilities.invokeLater(new ViewRunnable(viewName));
    }
  }

  public static View createView(String viewName) {

    if (viewName.equalsIgnoreCase(VN_VIRIAL)) {
      return new ApplicationView();
    }
    if (viewName.equalsIgnoreCase(VN_EMPTY_WIZARD)) {
      return new ClusterWizard(VN_EMPTY_WIZARD);
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
      v.initializeUI();
      v.display();
    }
  }
}