package simulate;

import java.awt.*;
import java.awt.event.*;
import java.beans.Beans;
import java.util.*;

public class Controller extends Container implements Runnable {

  public Integrator integrator;
  public Phase phase;
  public Simulation parentSimulation;
  Thread runner;
  private boolean initialized = false;
  private int maxSteps;

  public Controller() {
    setSize(100,40);
    maxSteps = Integer.MAX_VALUE;
  }

  public void add(Integrator i) {
//    super.add(i);
    this.integrator = i;
    i.parentController = this;

    for(Phase p=parentSimulation.firstPhase(); p!=null; p=p.nextPhase()) {
        i.registerPhase(p);
        p.gravity.addObserver(i);
        p.integrator = i;
    }
  }
  
  public void paint(Graphics g) {
//    initialized = true;    //assumes first call to paint comes when everything is in place
    if(Beans.isDesignTime()) {
      g.setColor(Color.red);
      g.drawRect(0,0,getSize().width-1,getSize().height-1);
      g.drawRect(1,1,getSize().width-3,getSize().height-3);
    }
    g.setColor(getBackground());
    paintComponents(g);
  }
  
//method to implement runnable
    public void start() {
        runner = new Thread(this);
        runner.start();
    }
    
    public int getMaxSteps() {return maxSteps;}
    public void setMaxSteps(int m) {maxSteps = m;}
    
    public void run() {
//        if(Beans.isDesignTime()) {return;}
//        while(!initialized) {;}   //stall until everything is in place
        integrator.initialize();
        integrator.run();
    }
}


