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
  
  public void registerPhases(Integrator i) {
    Enumeration e = phases.elements();
    while(e.hasMoreElements()) {
//      ((Phase)e.nextElement()).addPhaseIntegratorListener(pil);
    }
  }
  
  public void add(Display d) {
    super.add(d);
    if(haveIntegrator()) {
        controller.integrator.addIntegrationIntervalListener(d);
    }
    if(d.displayTool != null) {super.add(d.displayTool);}
    Phase p = (Phase)phases.firstElement();
    for(Meter m=p.firstMeter; m!=null; m=m.getNextMeter()) {d.setMeter(m);}
    d.repaint();
  }
  
  public void add(Phase p) {
    super.add(p);
    phases.addElement(p);
    p.parentSimulation = this;
    if(haveIntegrator()) {
        controller.integrator.registerPhase(p);
        p.gravity.addObserver(controller.integrator);
    }
  }
  
  public boolean haveIntegrator() {
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


