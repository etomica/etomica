package simulate;

import java.awt.*;
import java.awt.event.*;
import java.beans.Beans;

public class Controller extends Container implements Runnable {

  public Integrator integrator;
  public Phase phase;
  public Simulation parentSimulation;
  Thread runner;
  private boolean initialized = false;

  public Controller() {
    setSize(100,40);
  }

  public void add(Integrator i) {
    super.add(i);
    this.integrator = i;
    parentSimulation.registerPhases(i);
 //   runner = new Thread(this);
 //   runner.start();
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
    
    public void run() {
//        if(Beans.isDesignTime()) {return;}
//        while(!initialized) {;}   //stall until everything is in place
        integrator.initialize();
        integrator.run();
    }
}


