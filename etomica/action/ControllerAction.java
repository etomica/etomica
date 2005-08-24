package etomica.action;

import etomica.Controller;

 /**
  * Elementary action performed on a controller.
  */
public interface ControllerAction extends Action {

    public void setController(Controller c);
    
    public Controller getController();

}