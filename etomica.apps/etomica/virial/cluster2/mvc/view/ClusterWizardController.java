package etomica.virial.cluster2.mvc.view;

import javax.swing.SwingUtilities;

import etomica.virial.cluster2.mvc.Action;
import etomica.virial.cluster2.mvc.ActionResponse;
import etomica.virial.cluster2.mvc.ActionStatus;
import etomica.virial.cluster2.mvc.ViewResponse;
import etomica.virial.cluster2.mvc.ViewStatus;
import etomica.virial.cluster2.mvc.WizardController;
import etomica.virial.cluster2.mvc.WizardPageView;
import etomica.virial.cluster2.mvc.WizardView;

public class ClusterWizardController extends WizardController {

  @Override
  protected WizardView createWizard() {

    return new ClusterWizard("Cluster Creation Wizard");
  }

  @Override
  protected Action nextAction(ViewResponse response) {

    if (response.getStatus().isTerminated()) {
      return new DefaultWizardAction(ActionStatus.COMPLETE_SUCCESS, response.getView());
    }
    else {
      return new DefaultWizardAction(ActionStatus.CONTINUE_SUCCESS, response.getView());
    }
  }

  @Override
  protected WizardPageView nextPageView(ActionResponse actionResponse, ViewResponse viewResponse) {

    if (actionResponse == null || viewResponse == null) {
      return new ClusterWizardPage1();
    }
    // successful action
    if (actionResponse.getStatus() == ActionStatus.CONTINUE_SUCCESS) {
      // back button from the view
      if (viewResponse.getStatus() == ViewStatus.CONTINUE_PRIOR) {
        if (viewResponse.getView() instanceof ClusterWizardPage4) {
          return new ClusterWizardPage3();
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage3) {
          return new ClusterWizardPage2();
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage2) {
          return new ClusterWizardPage1();
        }
      }
      // next button from the view
      else if (viewResponse.getStatus() == ViewStatus.CONTINUE_NEXT) {
        if (viewResponse.getView() instanceof ClusterWizardPage3) {
          return new ClusterWizardPage4();
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage2) {
          return new ClusterWizardPage3();
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage1) {
          return new ClusterWizardPage2();
        }
      }
    }
    return new DefaultWizardPage();
  }

  public static void main(String[] args) {

    ApplicationUI.configure();
    SwingUtilities.invokeLater(new ClusterWizardController());
  }
}
