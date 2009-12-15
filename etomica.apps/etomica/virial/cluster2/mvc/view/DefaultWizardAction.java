package etomica.virial.cluster2.mvc.view;

import java.util.List;

import etomica.virial.cluster2.mvc.Action;
import etomica.virial.cluster2.mvc.ActionResponse;
import etomica.virial.cluster2.mvc.ActionStatus;
import etomica.virial.cluster2.mvc.MVCException;
import etomica.virial.cluster2.mvc.State;
import etomica.virial.cluster2.mvc.View;
import etomica.virial.cluster2.mvc.ViewResponse;

public class DefaultWizardAction implements Action {

  private ActionStatus status;
  private View view;

  public DefaultWizardAction(ActionStatus status, View view) {

    this.view = view;
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

      public View getView() {

        return view;
      }
    };
  }

}
