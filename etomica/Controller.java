package etomica;

import java.awt.Button;
import java.awt.Color;
import java.awt.Container;
import java.awt.Component;
import java.awt.event.*;
import java.util.*;

/**
 * Organizer of activities of the integrators.  Controller sets the protocol for
 * conduct of the simulation.  For example, the simulation might start and stop 
 * with the push of a button; it might run for a fixed number of cycles and exit;
 * it might run several simulations over a series of states.<p>
 * The Controller runs on its own thread, and it spawns Integrator processes
 * each on its own thread.
 */
public class Controller implements Simulation.Element, Runnable, java.io.Serializable, EtomicaElement {

  /**
   * List of integrators managed by the controller
   */
    protected final LinkedList integrators = new LinkedList();
  /**
   * Thread used to run the controller
   */
    transient Thread runner;
    private boolean initialized = false;
    private boolean autoStart = false;
    private int maxSteps;
    private String name;
    private Simulation parentSimulation;
    private boolean added = false;
    private Button startStopButton;
    
    private SimulationEventManager eventManager = new SimulationEventManager();

    public Controller() {
        this(Simulation.instance);
    }
    public Controller(Simulation sim) {
        parentSimulation = sim;
        maxSteps = Integer.MAX_VALUE;
        setMakeButton(true);
        parentSimulation.register(this);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Simple controller that enables start/stop of simulation using a button");
        return info;
    }

    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Controller.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
   /**
    * Resets the controller in some manner.
    * Performs no action for this controller, but may do something in subclasses
    */
    public void reset() {
        if(startStopButton != null) startStopButton.reset();
    }
 
   /**
    * Adds an integrator to the list of integrators managed by this controller.
    * Starts the integrator immediately if autoStart is true
    */
    public void add(Integrator i) {
        integrators.add(i);
        i.setController(this);
        i.setMaxSteps(maxSteps);
        if(autoStart) i.start();
    }
    
    /**
     * @return a list of the integrators managed by this controller
     */
    public LinkedList integrators() {return integrators;}
    
    /**
     * Accessor method for the autoStart flag.
     * autoStart flags for whether the integrators start as soon as added to the controller.
     * Default is false
     */
    public void setAutoStart(boolean b) {autoStart = b;}
    /**
     * Accessor method for the autoStart flag.
     * autoStart flags for whether the integrators start as soon as added to the controller.
     * Default is false
     */
    public boolean getAutoStart() {return autoStart;}
    
    /**
     * Accessor for maxSteps, the total number of steps performed by each integrator.
     */
    public int getMaxSteps() {return maxSteps;}
    /**
     * Accessor for maxSteps, the total number of steps performed by each integrator.
     */
    public void setMaxSteps(int m) {
        maxSteps = m;
        for(java.util.Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            ((Integrator)iter.next()).setMaxSteps(m);
        }
    }

    /**
     * Method to start this controller's thread, causing execution of the run() method
     */
    public void start() {
        if(runner != null) return;
        parentSimulation().elementCoordinator.go();
        runner = new Thread(this);
        runner.start();
    }
                    
    /**
     * Activity performed by this controller once its start method is called.
     * This controller simply initiates each integrator on its own thread, which
     * is initiated by calling its start method.  Controller subclasses would override
     * this method to perform some other function.<p>
     * An integrator can be run on the controller's thread by instead calling
     * the run() method of the integrator.
     */
    public void run() {
        parentSimulation().elementCoordinator.go();  //perhaps redundant, but safe since it returns without performing any action if already completed
        for(java.util.Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            ((Integrator)iter.next()).start();
        }
        runner = null;
    }
    
    /**
     * Sends a halt (stop) signal to all integrators.
     */
    public void halt() {
        for(java.util.Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            ((Integrator)iter.next()).halt();
        }
    }        
    
    public synchronized void addListener(SimulationEventListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(SimulationEventListener listener) {
        eventManager.removeListener(listener);
    }

    protected void fireEvent(SimulationEvent event) {
        eventManager.fireEvent(event);
    }    
    
    private boolean makeButton = false;
    
    public void setMakeButton(boolean b) {
        if(!b || makeButton) return; //return if argument is false or already made button
        startStopButton = new Button();
        makeButton = true;
    }
    public boolean getMakeButton() {return makeButton;}
    
    public Button getButton() {return startStopButton;}
    
    /**
     * Accessor method of the name of this object
     * 
     * @return The given name
     */
    public final String getName() {return name;}

    /**
     * Method to set the name of this object
     * 
     * @param name The name string to be associated with this object
     */
    public final void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the object
     */
    public String toString() {return getName();}  //override Object method
          

    /**
     * A device that presents a button to pause and resume all integrators
     * operated by this controller.  This is a member class of this controller.
     * An instance of it is created by the syntax<br>
     *   Controller.Button button = controller.new Button()<br>
     * where "controller" is the name of the instance of this controller
     */
    public class Button extends Device implements ActionListener {
        
        /**
         * The awt Button that the user interacts with to control the Controller
         */
        public javax.swing.JButton button;
        private boolean firstClick = true;
        private boolean running = false;

        public Button() {
            super(Controller.this.parentSimulation());
            if(button == null) button = new javax.swing.JButton("Start");
 //           button.setSize(60,40);
            autoStart = false;
            button.addActionListener(this);  //register this as a listener to the awt Button
        }
        
        /**
         * Sets the button to the state when it is first created, and before it has
         * started the controller
         */
        public void reset() {
            firstClick = true;
            button.setText("Start");
        }
        /**
         * A method that mimics the action of a button click.  Useful if another
         * object wants to pause the simulation in a way that leaves the button in
         * a consistent state
         */
        public void clickButton() {  
            actionPerformed(new ActionEvent(button,0,""));
        }
        
        /**
         * A method that mimics the action of a button click, but only if the integrators are not already paused.
         * If the integrators are currently paused, the method has no action.
         */
        public void clickForPause() {  //clicks button to pause if not already paused
            if(running) clickButton();
        }
              
        /**
         * A method that mimics the action of a button click, but only if the integrators are already paused.
         * If the integrators are not currently paused, the method has no action.
         */
        public void clickForUnpause() {  //clicks button to unpause if already paused
            if(!running) clickButton();
        }

        public void setButtonColor(Color c) {button.setBackground(c);}
        public Color getButtonColor() {return button.getBackground();}
              
        public void setButtonTextColor(Color c) {button.setForeground(c);}
        public Color getButtonTextColor() {return button.getForeground();}
        
        /**
         * ActionListener interface method (button is a self-listener).  
         * First click invokes the controller's start() method, while subsequent clicks 
         * toggle between calling pause() and unpause() methods of all integrators.
         * Label on button is modified accordingly.
         */
        public void actionPerformed(ActionEvent evt) {
            if(firstClick) {
                start();
                running = true;
                firstClick = false;
                button.setText("Pause");
//                etomica.gui.EtomicaToolBar.startButton.setEnabled(false);
            }
            else if(running) {
                running = false;
                for(java.util.Iterator iter=integrators.iterator(); iter.hasNext(); ) {
                    ((Integrator)iter.next()).pause();
                }
                button.setText("Continue");
            }
            else {
                for(java.util.Iterator iter=integrators.iterator(); iter.hasNext(); ) {
                    ((Integrator)iter.next()).unPause();
                }
                running = true;
	            button.setText("Pause");
	        }
	    }
        
        /**
         * Simulation.GraphicElement method that returns a handle to the awt Button instance used by this class.
         */
        public Component graphic(Object obj) {
            if(button == null) button = new javax.swing.JButton("Start");
            return button;
        }
            
    }//end of Controller.Button
}


