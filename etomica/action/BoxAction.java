package etomica.action;

import etomica.box.Box;

 /**
  * Elementary action performed on a box.
  */
public interface BoxAction extends Action {

    public void setBox(Box box);
    
    public Box getBox();

}