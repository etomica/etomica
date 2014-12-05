/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import static etomica.virial.cluster2.mvc.view.ClusterWizardState.DEFVAL_MONOCHROMATIC;
import static etomica.virial.cluster2.mvc.view.ClusterWizardState.KEY_COLOR_SCHEME;
import static etomica.virial.cluster2.mvc.view.ClusterWizardState.KEY_ISOMORPH_FREE;
import etomica.virial.cluster2.mvc.*;

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
  protected WizardPageView createPageView(ActionResponse actionResponse, ViewResponse viewResponse) {

    if (actionResponse == null || viewResponse == null) {
      return new ClusterWizardPage1(this);
    }
    // successful action
    if (actionResponse.getStatus() == ActionStatus.CONTINUE_SUCCESS) {
      // back button from the view
      if (viewResponse.getStatus() == ViewStatus.CONTINUE_PRIOR) {
        if (viewResponse.getView() instanceof ClusterWizardPage5) {
          // do we show the color assignment page?
          if ((Boolean) getState().getProperty(KEY_COLOR_SCHEME).equals(DEFVAL_MONOCHROMATIC)) {
            return new ClusterWizardPage3(this);
          }
          else {
            return new ClusterWizardPage4(this);
          }
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage4) {
          return new ClusterWizardPage3(this);
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage3) {
          return new ClusterWizardPage2(this);
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage2) {
          return new ClusterWizardPage1(this);
        }
      }
      // next button from the view
      else if (viewResponse.getStatus() == ViewStatus.CONTINUE_NEXT) {
        if (viewResponse.getView() instanceof ClusterWizardPage4) {
          return new ClusterWizardPage5(this);
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage3) {
          // do we show the color assignment page?
          if ((Boolean) getState().getProperty(KEY_COLOR_SCHEME).equals(DEFVAL_MONOCHROMATIC)) {
            return new ClusterWizardPage5(this);
          }
          else {
            return new ClusterWizardPage4(this);
          }
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage2) {
          return new ClusterWizardPage3(this);
        }
        else if (viewResponse.getView() instanceof ClusterWizardPage1) {
          return new ClusterWizardPage2(this);
        }
      }
    }
    return new ClusterWizardPage1(this);
  }

  @Override
  protected WizardState createWizardState() {

    return new ClusterWizardState();
  }
}
