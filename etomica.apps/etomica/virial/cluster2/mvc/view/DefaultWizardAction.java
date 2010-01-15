package etomica.virial.cluster2.mvc.view;

import java.util.List;

import etomica.virial.cluster2.mvc.*;

public class DefaultWizardAction implements Action {

  private ActionStatus status;

  public DefaultWizardAction(ActionStatus status) {

    this.status = status;
  }

  public ActionResponse execute(final ViewResponse response) {

    final Action thisAction = this;

    return new ActionResponse() {

      public Action getAction() {

        return thisAction;
      }

      public ActionStatus getStatus() {

        return status;
      }

      public List<MVCException> getErrors() {

        return null;
      }
    };
  }

}
