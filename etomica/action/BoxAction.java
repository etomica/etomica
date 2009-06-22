package etomica.action;

import etomica.api.IBox;

 /**
  * Elementary action performed on a box.
  */
public interface BoxAction extends IAction {

    public void setBox(IBox box);
    
    public IBox getBox();

}