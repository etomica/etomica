package simulate;

import java.awt.*;
import java.awt.event.*;
import java.beans.Beans;
import java.util.*;

public class Simulation extends Container {

  public static int D;  //dimension (2-D, 3-D, etc;)

  public Controller controller;
  PhaseSpace firstPhaseSpace;
  PhaseSpace lastPhaseSpace;
  Display firstDisplay;
  Display lastDisplay;
  
  public Simulation() {
    this(2);
  }
  
  public Simulation(int D) {
    setD(D);
    setSize(400,300);
  }
  
  //For property sheet (D is static)
  public void setD(int d) {D = d;}
  public int getD() {return D;}


  public void add(Controller c) {
    if(controller != null) {return;}  //already added a controller
    super.add(c);
    controller = c;
    c.parentSimulation = this;
  }
    
  public void add(Display d) {
    super.add(d);
    d.parentSimulation = this;
    if(lastDisplay != null) {lastDisplay.setNextDisplay(d);}
    else {firstDisplay = d;}
    lastDisplay = d;
    if(haveIntegrator()) {
        controller.integrator.addIntegrationIntervalListener(d);
    }
    if(d.displayTool != null) {super.add(d.displayTool);}
    for(PhaseSpace p=firstPhaseSpace; p!=null; p=p.nextPhaseSpace()) {
        d.setPhaseSpace(p);
    }
    d.repaint();
  }
  
  public void add(PhaseSpace p) {
    super.add(p);
    p.parentSimulation = this;
    if(lastPhaseSpace != null) {lastPhaseSpace.setNextPhase(p);}
    else {firstPhaseSpace = p;}
    lastPhaseSpace = p;
    if(haveIntegrator()) {
        controller.integrator.registerPhaseSpace(p);
        p.gravity.addObserver(controller.integrator);
        p.integrator = controller.integrator;
    }
    for(Display d=firstDisplay; d!=null; d=d.getNextDisplay()) {
        d.setPhaseSpace(p);
    }
  }
  
  public PhaseSpace firstPhaseSpace() {return firstPhaseSpace;}
  public PhaseSpace lastPhaseSpace() {return lastPhaseSpace;}
  
  public boolean haveIntegrator() {
    return (controller != null && controller.integrator != null);
  }
  
/*  public void paint(Graphics g) {
    if(Beans.isDesignTime()) {
      g.setColor(Color.red);
      g.drawRect(0,0,getSize().width-1,getSize().height-1);
      g.drawRect(1,1,getSize().width-3,getSize().height-3);
    } 
    g.setColor(getBackground());
    paintComponents(g);
  }
  */
}


