package simulate;
import java.awt.*;
import java.beans.Beans;
import java.awt.event.*;

    public abstract class Display extends Panel implements simulate.IntegrationIntervalListener {

	Simulation parentSimulation;
	Phase phase;
    int updateInterval;
    int iieCount;
//    Component displayTool = null;  //displayTool is some component inside a Display object that is doing all the painting (not often used)
    private Display nextDisplay;
    private Display previousDisplay;
        
    public void addNotify() {
        super.addNotify();
        Container p = getParent();
        if(p instanceof Simulation) {return;}
        while(p != null) {
            p = p.getParent();
            if(p instanceof Simulation) {
                ((Simulation)p).addDisplay(this);
                System.out.println("OK in display");
                return;
            }
        }
        System.out.println("Warning:  Display has no parent Simulation");
    }            
    
    public Display () {
	    setUpdateInterval(1);
    }
    
    public final Display nextDisplay() {return nextDisplay;}
    public final Display getPreviousDisplay() {return previousDisplay;}
   /**
    * Sets the display following this one in the linked list of displays.
    *
    * @param d the display to be designated as this display's nextDisplay
    */
    public final void setNextDisplay(Display d) {
      this.nextDisplay = d;
      d.previousDisplay = this;
    }
    
    public void setPhase(Phase p) {phase = p;}
    
    public abstract void doUpdate();
        
    public void integrationIntervalAction(IntegrationIntervalEvent evt) {
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        doUpdate();
            repaint();
	    }
    }

	public void setParentSimulation(Simulation s) {parentSimulation = s;}
	public Simulation parentSimulation() {return parentSimulation;}

    public final int getUpdateInterval() {return updateInterval;}
    public final void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
    }
}
