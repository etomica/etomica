package etomica.action;

import etomica.api.IBox;
import etomica.box.Box;

 /**
  * Elementary action performed on a box.
  */
public interface BoxAction extends Action {

    public void setBox(Box box);
    
    public IBox getBox();

}