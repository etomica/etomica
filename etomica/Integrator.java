package simulate;
import java.io.*;
import java.awt.event.*;
import java.awt.*;
import java.util.*;

public abstract class Integrator extends Container implements PhaseIntegratorListener, Serializable, MouseListener, Runnable {

    //VARIABLES used to make runnable    
  Thread runner;
  int running = 0;
 
  Phase firstPhase;
  Phase[] phase;
  int nPhases = 0;
  int nPhasesMax = 1;
  protected int sleepPeriod = 10;
  protected transient Vector listeners=null;
  private Vector integrationIntervalListeners;
  int integrationInterval;  // number of drawTimeSteps between IntegrationIntervalEvent firing
  int integrationCount;
  double drawTimeStep;
  private int neighborListUpdateInterval;

  public double temperature = 300;
  public boolean isothermal = false;

  public Integrator() {
    integrationIntervalListeners = new Vector();
    integrationInterval = 10;
    integrationCount = 0;
    drawTimeStep = 0.0005;
    neighborListUpdateInterval = Integer.MAX_VALUE;
    phase = new Phase[nPhasesMax];
  }
  
  public final int getSleepPeriod() {return sleepPeriod;}
  public final void setSleepPeriod(int s) {sleepPeriod = s;}

// abstract methods
  public abstract void doStep(double tStep);
  public abstract void initialize();
  
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
  
  public final int getNeighborListUpdateInterval() {return neighborListUpdateInterval;}
  public final void setNeighborListUpdateInterval(int interval) {neighborListUpdateInterval = interval;}
  
  // Event-related methods
  public void phaseIntegratorNotify(PhaseIntegratorEvent pie) {
    if(nPhases == nPhasesMax) {return;}
    phase[nPhases] = (Phase)pie.getSource();
    nPhases++;
    firstPhase = phase[0];
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
        listener.updateAverage(iie);
    }
  }
 
//methods to implement runnable    
    public void start() {;}

    public void stop() {
      if(runner!=null) {
          runner.stop();
          runner = null;
      }
    }

    public void run() {
        int ic = 0;
        while(true) {
            this.doStep(drawTimeStep);
            integrationCount++;
            if(integrationCount == Integer.MAX_VALUE) {integrationCount = 1;}
            if(integrationCount % neighborListUpdateInterval == 0) {
                firstPhase.updateNeighbors();
            }
            if(integrationCount % integrationInterval == 0) {
                fireIntegrationIntervalEvent(new IntegrationIntervalEvent(this,firstPhase));
            }
            for (int i=0; i<nPhases; i++) {phase[i].repaint();}

            try { Thread.sleep(sleepPeriod); }
            catch (InterruptedException e) { }
        }
    }

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
    
}

