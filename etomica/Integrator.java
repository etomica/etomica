package simulate;
import java.io.*;
import java.awt.event.*;
import java.awt.*;
import java.util.*;

public abstract class Integrator extends Container implements Observer, Serializable, MouseListener, Runnable {

    //VARIABLES used to make runnable    
  Thread runner;
  int running = 0;
  protected int maxSteps = Integer.MAX_VALUE;
 
  Phase firstPhase;
  Phase[] phase;
  int phaseCount = 0;
  int phaseCountMax = 1;
  protected int sleepPeriod = 10;
  protected transient Vector listeners=null;
  private Vector integrationIntervalListeners;
  int integrationInterval;  // number of drawTimeSteps between IntegrationIntervalEvent firing
  int integrationCount;
  double drawTimeStep;
  boolean doSleep = true;
  private int neighborListUpdateInterval;

  public double temperature = 300;
  public boolean isothermal = false;

  public Integrator() {
    integrationIntervalListeners = new Vector();
    integrationInterval = 10;
    integrationCount = 0;
    drawTimeStep = 0.0005;
    neighborListUpdateInterval = Integer.MAX_VALUE;
    phase = new Phase[phaseCountMax];
  }
  
  public abstract Agent makeAgent(Atom a);
  
  protected void deployAgents() {  //puts an Agent of this integrator in each atom of all phases
    for(Phase p=firstPhase; p!=null; p=p.nextPhase()) {
        for(Atom a=p.firstAtom(); a!=null; a=a.nextAtom()) {
            a.setIntegratorAgent(makeAgent(a));
        }
    }
  }
    
  public final int getSleepPeriod() {return sleepPeriod;}
  public final void setSleepPeriod(int s) {sleepPeriod = s;}

// abstract methods
  public abstract void doStep(double tStep);
  public abstract void initialize();  //put here a call to deployAgents, among other things
  
  // Introspected properties
  public final void setTemperature(double t) {temperature = t;}
  public final void setITemperature(int t) {temperature = (double)t;}
  public final double getTemperature() {return temperature;}

  public void setIsothermal(boolean b) {isothermal = b;}
  public boolean isIsothermal() {return isothermal;}
    
  public final double getDrawTimeStep() {return drawTimeStep;}
  public final void setDrawTimeStep(double t) {drawTimeStep = t;}
  
  public final int getIntegrationInterval() {return integrationInterval;}
  public final void setIntegrationInterval(int interval) {integrationInterval = interval;}
  
//  public final int getNeighborListUpdateInterval() {return neighborListUpdateInterval;}
//  public final void setNeighborListUpdateInterval(int interval) {neighborListUpdateInterval = interval;}
    
  public void registerPhase(Phase p) {
    if(phaseCount == phaseCountMax) {return;}
    phase[phaseCount] = p;
    phaseCount++;
    firstPhase = phase[0];
    for(Meter m=p.firstMeter; m!=null; m=m.nextMeter()) {
        addIntegrationIntervalListener(m);
    }
  }
  
  public synchronized void addIntegrationIntervalListener(IntegrationIntervalListener iil) {
    integrationIntervalListeners.addElement(iil);
  }

  public synchronized void removeIntegrationIntervalListener(IntegrationIntervalListener iil) {
    integrationIntervalListeners.removeElement(iil);
  }

  public void fireIntegrationIntervalEvent(IntegrationIntervalEvent iie) {
    Vector currentListeners = null;
    synchronized(this){
        currentListeners = (Vector)integrationIntervalListeners.clone();
    }
    for(int i = 0; i < currentListeners.size(); i++) {
        IntegrationIntervalListener listener = (IntegrationIntervalListener)currentListeners.elementAt(i);
        listener.integrationIntervalAction(iie);
    }
  }
  
  /**
   Update method for Observer interface
   */
  public void update(Observable o, Object arg) {}
 
//methods to implement runnable    
    public void start() {;}

    public void stop() {
      if(runner!=null) {
          runner.stop();
          runner = null;
      }
    }

    public int getMaxSteps() {return maxSteps;}
    public void setMaxSteps(int m) {maxSteps = m;}
    
    public void run() {
        int ic = 0;
        int nSteps = 0;
        while(nSteps < maxSteps) {
            this.doStep(drawTimeStep);
            integrationCount++;
            if(integrationCount == Integer.MAX_VALUE) {integrationCount = 1;}
///            if(integrationCount % neighborListUpdateInterval == 0) {
///                firstPhase.updateNeighbors();
///            }
            if(integrationCount % integrationInterval == 0) {
//                for(Phase p=firstPhase; p!=null; p=p.nextPhase()) {p.space.repositionMolecules();}                
                fireIntegrationIntervalEvent(new IntegrationIntervalEvent(this,firstPhase));
            }
//            for (int i=0; i<phaseCount; i++) {phase[i].repaint();}  //not needed now that displayConfiguration does painting
            if(doSleep) {
                try { Thread.sleep(sleepPeriod); }
                catch (InterruptedException e) { }
            }
            nSteps++;
//            System.out.println(nSteps);
        }
        System.exit(0);
    }
    
    public void setDoSleep(boolean b) {doSleep = b;}
    public boolean isDoSleep() {return doSleep;}

//methods to implement mouseListener
    public void mouseClicked(MouseEvent evt) {
        if(evt.getSource() instanceof Button) {
            Button button = (Button)evt.getSource();
            if(runner==null) {
                this.initialize();
                runner = new Thread(this);
                runner.start();
                running = 1;
   //             button.setBackground(Color.red);
                button.setLabel("Pause");
            }
            else if(running==1) {
                runner.suspend();
                running = 0;
   //             button.setBackground(Color.green);
                button.setLabel("Continue");
            }
            else {
                runner.resume();
                running = 1;
    //	        button.setBackground(Color.red);
	            button.setLabel("Pause");
	        }
	    }
	}
    
    public void mouseEntered(MouseEvent e) {
 //       mouseClicked(e);
    }    
    public void mouseExited(MouseEvent e) {
 //       mouseClicked(e);
    }    
    public void mousePressed(MouseEvent e) {
  //      mouseClicked(e);
    }    
    public void mouseReleased(MouseEvent e) {
 //       mouseClicked(e);
    }
 
// Class generated by integrator as one of the properties of each atom

    interface Agent {}
    
}

