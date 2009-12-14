package etomica.virial.cluster2.mvc.view;

import java.util.List;

import etomica.virial.cluster2.mvc.Action;
import etomica.virial.cluster2.mvc.ActionResponse;
import etomica.virial.cluster2.mvc.ActionStatus;
import etomica.virial.cluster2.mvc.MVCException;
import etomica.virial.cluster2.mvc.State;
import etomica.virial.cluster2.mvc.ViewResponse;

public class DefaultWizardAction implements Action {

  private ActionStatus status;

  public DefaultWizardAction(ActionStatus status) {

    this.status = status;
  }

  public ActionResponse execute(final ViewResponse response, final State state) {

    final Action thisAction = this;

    return new ActionResponse() {

      public Action getAction() {

        return thisAction;
      }

      public ActionStatus getStatus() {

        return status;
      }

      public State getData() {

        return state;
      }

      public List<MVCException> getErrors() {

        return null;
      }
    };
  }

}
