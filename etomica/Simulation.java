package simulate;

import java.awt.*;
import java.awt.event.*;
import java.beans.Beans;
import java.util.*;

public class Simulation extends Container {

  public Controller controller;
  public Vector phases = new Vector();
  
  public Simulation() {
    setSize(400,300);
  }

  public void add(Controller c) {
    if(controller != null) {return;}  //already added a controller
    super.add(c);
    controller = c;
    c.parentSimulation = this;
  }
  
  public void registerPhases(PhaseIntegratorListener pil) {
    Enumeration e = phases.elements();
    while(e.hasMoreElements()) {
      ((Phase)e.nextElement()).addPhaseIntegratorListener(pil);
    }
  }
  
  public void add(Phase p) {
    super.add(p);
    phases.addElement(p);
    if(haveIntegrator()) {
        p.addPhaseIntegratorListener(controller.integrator);
        p.gravity.addObserver(controller.integrator);
    }
  }
  
  private boolean haveIntegrator() {
    return (controller != null && controller.integrator != null);
  }
  
  public void paint(Graphics g) {
    if(Beans.isDesignTime()) {
      g.setColor(Color.red);
      g.drawRect(0,0,getSize().width-1,getSize().height-1);
      g.drawRect(1,1,getSize().width-3,getSize().height-3);
    } 
    g.setColor(getBackground());
    paintComponents(g);
  }
}


