package etomica.action;

import etomica.action.activity.IController;
import etomica.api.IAction;

 /**
  * Elementary action performed on a controller.
  */
public interface ControllerAction extends IAction {

    public void setController(IController c);
    
    public IController getController();

}