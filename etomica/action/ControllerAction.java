package etomica.action;

import etomica.api.IAction;
import etomica.api.IController;

 /**
  * Elementary action performed on a controller.
  */
public interface ControllerAction extends IAction {

    public void setController(IController c);
    
    public IController getController();

}