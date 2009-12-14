package etomica.virial.cluster2.mvc.view;

import javax.swing.SwingUtilities;

import etomica.virial.cluster2.mvc.Action;
import etomica.virial.cluster2.mvc.ActionResponse;
import etomica.virial.cluster2.mvc.ActionStatus;
import etomica.virial.cluster2.mvc.ViewResponse;
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
      return new DefaultWizardAction(ActionStatus.COMPLETE_SUCCESS);
    }
    else {
      return new DefaultWizardAction(ActionStatus.CONTINUE_SUCCESS);
    }
  }

  @Override
  protected WizardPageView nextPageView(ActionResponse response) {

    return new DefaultWizardPage();
  }

  public static void main(String[] args) {

    ApplicationUI.configure();
    SwingUtilities.invokeLater(new ClusterWizardController());
  }
}
