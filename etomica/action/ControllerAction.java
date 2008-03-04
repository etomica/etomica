package etomica.action;

import etomica.api.IController;

 /**
  * Elementary action performed on a controller.
  */
public interface ControllerAction extends Action {

    public void setController(IController c);
    
    public IController getController();

}