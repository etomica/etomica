//This class includes a main method to demonstrate its use
package etomica;

import javax.swing.*;
//import java.util.*;
import etomica.units.UnitSystem;
import java.beans.Beans;

//Java2 imports
//import java.util.HashMap;
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.HashMap;
import etomica.utility.LinkedList;
import etomica.utility.Iterator;

/**
 * The main class that organizes the elements of a molecular simulation.
 * Holds a single space object that is referenced in
 * many places to obtain spatial elements such as vectors and boundaries.  Also
 * holds an object that specifies the unit system used to default all I/O.  A single
 * instance of Simulation is held as a static field, and which forms the default
 * Simulation class needed in the constructor of all simulation elements.
 *
 * @author David Kofke
 */
public class Simulation extends javax.swing.JPanel implements java.io.Serializable {
    
    public String getVersion() {return "Simulation:01.06.05";}
    /**
     * Flag indicating whether simulation is being run within Etomica editor application.
     * This is set to true by Etomica if it is running; otherwise it is false.
     */
    public static boolean inEtomica = false;
    /**
     * Class that implements the final tying up of the simulation elements before starting the simulation.
     * Default choice is the CoordinatorOneIntegrator.
     */
    public Mediator elementCoordinator;
    private HashMap elementLists = new HashMap(16);
    private final PotentialMaster potentialMaster = new PotentialMaster(this);
    
   /**
    * Object describing the nature of the physical space in which the simulation is performed
    */
    public final Space space; //would like to make final, but compiler doesn't allow
    
    /**
     * List of all controllers that have been instantiated.
     */
    public LinkedList controllerList = new LinkedList();
    /**
     * List of all phases that have been instantiated.
     */
    public LinkedList phaseList = new LinkedList();
    /**
     * List of all species that have been instantiated.
     */
    public LinkedList speciesList = new LinkedList();
    /**
     * List of all displays that have been instantiated.
     */
    public LinkedList displayList = new LinkedList();
    /**
     * List of all devices that have been instantiated.
     */
    public LinkedList deviceList = new LinkedList();
    /**
     * List of all integrators that have been instantiated.
     */
    public LinkedList integratorList = new LinkedList();
    /**
     * List of all potentials that have been instantiated.
     */
    public LinkedList potentialList = new LinkedList();
    /**
     * List of all meters that have been instantiated.
     */
    public LinkedList meterList = new LinkedList();
    /**
     * List of all simulation elements.
     */
     private LinkedList allElements = new LinkedList();

    //default unit system for I/O (internal calculations are all done in simulation units)
    private static UnitSystem unitSystem = new UnitSystem.Sim();
    
	public final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
	public final javax.swing.JPanel displayBoxPanel = new JPanel(new java.awt.GridBagLayout());
//    public final javax.swing.JPanel devicePanel = new JPanel(new java.awt.GridLayout(0,1),false);
    public final javax.swing.JPanel devicePanel = new JPanel(new java.awt.GridBagLayout());
//    public final javax.swing.JPanel devicePanel = new JPanel(new java.awt.FlowLayout());

    /**
     * A static instance of a Simulation, for which the current value at any time is
     * used as a default simulation in many places.
     */
    public static Simulation instance = new Simulation(new Space2D());
       
    public Simulation() {
        this(new Space2D());
    }
    
    /**
     * Constructor requires specification of the space used by the simulation
     */
    public Simulation(Space s) {
        super();
        space = s;
        setName("Simulation" + Integer.toString(instanceCount++));
        elementLists.put(Potential.class, potentialList);
        elementLists.put(Species.class, speciesList);
        elementLists.put(Integrator.class, integratorList);
        elementLists.put(Phase.class, phaseList);
        elementLists.put(Controller.class, controllerList);
        elementLists.put(Display.class, displayList);
        elementLists.put(MeterAbstract.class, meterList);
        elementLists.put(Device.class, deviceList);
        elementCoordinator = new Mediator(this);
        setSize(800,550);
        setLayout(new java.awt.FlowLayout());
        add(devicePanel);
        add(displayBoxPanel);
        add(displayPanel);
/*        setLayout(new java.awt.BorderLayout());
        add(devicePanel, java.awt.BorderLayout.NORTH);
        add(displayPanel, java.awt.BorderLayout.EAST);
        add(displayBoxPanel, java.awt.BorderLayout.WEST);*/
/*        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc = new java.awt.GridBagConstraints();
        gbc.gridx = 0;
        add(devicePanel, gbc);
        add(displayBoxPanel, gbc);
        gbc.gridx = 1;
        add(displayPanel, gbc);*/
        //workaround for JTabbedPane bug in JDK 1.2
        displayPanel.addChangeListener(
            new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent event) {
                    displayPanel.invalidate();
                    displayPanel.validate();
                }
        });
    }//end of constructor
    
    /**
     * Accessor method for the default I/O unit system.
     */
    public static final UnitSystem unitSystem() {return unitSystem;}
    /**
     * Accessor method for the default I/O unit system.
     */
    public static final void setUnitSystem(UnitSystem us) {unitSystem = us;}
              
    public final Space space() {return space;}
    /**
     * @return the <code>nth</code> instantiated phase (indexing from zero)
     */
    public final Phase phase(int n) {return (Phase)phaseList.get(n);}
    /**
     * @return the <code>nth</code> instantiated species (indexing from zero)
     */
    public final Species species(int n) {return (Species)speciesList.get(n);}
    /**
     * @return the <code>nth</code> instantiated potential (indexing from zero)
     */
    public final Potential potential(int n) {return (Potential)potentialList.get(n);}
    /**
     * @return the <code>nth</code> instantiated controller (indexing from zero)
     */
    public final Controller controller(int n) {return (Controller)controllerList.get(n);}
    /**
     * @return the <code>nth</code> instantiated integrator (indexing from zero)
     */
    public final Integrator integrator(int n) {return (Integrator)integratorList.get(n);}
    /**
     * @return the <code>nth</code> instantiated meter (indexing from zero)
     */
    public final MeterAbstract meter(int n) {return (MeterAbstract)meterList.get(n);}
    /**
     * @return the <code>nth</code> instantiated display (indexing from zero)
     */
    public final Display display(int n) {return (Display)displayList.get(n);}
    /**
     * @return the <code>nth</code> instantiated device (indexing from zero)
     */
    public final Device device(int n) {return (Device)deviceList.get(n);}
    
    public final LinkedList phaseList() {return phaseList;}
    public final LinkedList meterList() {return meterList;}
    public final LinkedList speciesList() {return speciesList;}
    public final LinkedList integratorList() {return integratorList;}
    public final LinkedList controllerList() {return controllerList;}
    public final LinkedList potentialList() {return potentialList;}
    public final LinkedList displayList() {return displayList;}
    public final LinkedList deviceList() {return deviceList;}
  
    public void register(Simulation.Element element) {
        if(element == potentialMaster || potentialMaster == null) return;
        LinkedList list = (LinkedList)elementLists.get(element.baseClass());
        if(list.contains(element)) return;
        if(element instanceof P2IdealGas) return;
        element.setName(element.getClass().getName().substring(8) + Integer.toString(list.size()));
        list.add(element);
        allElements.add(element);
    }//end of register method
                
    public void unregister(Simulation.Element element) {
        LinkedList list = (LinkedList)elementLists.get(element.baseClass());
        if(!list.contains(element)) return;
        list.remove(element);
        allElements.remove(element);
    }
                
    /**
     * Method invoked in the constructor of a Display object to list it with the simulation
     */
    public void unregister(Display d) {
        if(!displayList.contains(d)) return;
        displayList.remove(d);
        allElements.remove(d);
        integrator(0).removeIntervalListener(d);
        for (int i = 0; i < getComponentCount(); i++) {
	        if (getComponent(i).getName() == "displayPanel"){
	            javax.swing.JTabbedPane tp = ((javax.swing.JTabbedPane)getComponent(i));
	            tp.remove(tp.getSelectedComponent());
	            tp.repaint();
	            elementCoordinator.completed = false;
	            
	        }
	    }
    }
    
    /**
     * Method invoked in the constructor of a Device object to list it with the simulation
     */
    public void unregister(Device d) {
        if(!deviceList.contains(d)) return;
        deviceList.remove(d);
        allElements.remove(d);
        for (int i = 0; i < getComponentCount(); i++) {
	        if (Simulation.instance.getComponent(i).getName() == "devicePanel"){
                javax.swing.JPanel dp = ((javax.swing.JPanel)getComponent(i));
                dp.remove(dp.getComponent(1));
                dp.repaint();
                elementCoordinator.completed = false;
            }
        }
    }
              
                  
    /**
     * Method invoked in the constructor of a Species object to list it with the simulation
     */
    public void unregister(Species species) {
        if(!speciesList.contains(species)) return;
        //remove from manager
        speciesList.remove(species);
        allElements.remove(species);
    }
                
                
    /**
     * Method invoked in the constructor of a Potential2 object to list it with the simulation
     */
/*    public void unregister(Potential2 p2) {
        if(!potential2List.contains(p2)) return;
        if(p2 instanceof P2IdealGas) return;
        potential2List.remove(p2);
        allElements.remove(p2);
    }
*/ 
    public LinkedList allElements() {
        LinkedList list;
 //       synchronized(this) {
            list = (LinkedList)allElements.clone();
 //       }
        return list;
    }
    
    public final PotentialGroup potentialMaster() {return potentialMaster;}
            
    /**
     * Accessor of the first species in the linked list of species
     */
    public final Species firstSpecies() {return firstSpecies;}
    /**
     * Accessor of the last species in the linked list of species
     */
    public final Species lastSpecies() {return lastSpecies;}
    /**
     * @return the number of species that have been constructed and registered with the simulation
     */
     public final int speciesCount() {return speciesList.size();}
    
    private static int instanceCount = 0;
    private String name;
    public void setName(String newName) {name = newName;}
    public String getName() {return name;}
    public String toString() {return getName();}
    
    /**
    * Symmetric array of all two-body potentials.  Potentials are associated with species, and each species
    * is assigned a unique index to idenfity it.  Potential2[i][j] is the two-body potential
    * for Species indexed i and j, respectively.  The potential for i=j is merely the one describing the 
    * molecular interactions for that species.
    * 
    * @see Species#speciesIndex
    * @see Potential2
    */
    public Potential2[][] potential2;
      
    /**
    * Array of all one-body potentials.  Potentials are associated with species, and each species
    * is assigned a unique index to idenfity it.  Potential1[i] is the one-body potential
    * for Species indexed i.
    * 
    * @see Species#speciesIndex
    * @see Potential1
    */
    public Potential1[] potential1;
      
    /**
    * First species in the linked list of species in this phase.
    */
    public Species firstSpecies;
     
    /**
    * Last species in the linked list of species in this phase.
    */
    public Species lastSpecies;
    
    public static final java.util.Random random = new java.util.Random();
        
    /**
     * A visual display of the simulation via a JPanel.
     * This may become more important if Simulation itself is revised to not extend JPanel.
     */
     public javax.swing.JPanel panel() {return this;}
     
     /**
      * Returns the mediator that coordinates the elements of the simulation.
      * The is the same as the elementCoordinator field, but provides another
      * way to access it.  This method may someday supercede direct access to
      * the elementCoordinator field, so it is the preferred way to access it.
      */
     public Mediator mediator() {return elementCoordinator;}
     
   /**
    * A marker interface to indicate that the class is an element of a Simulation
    */
    public interface Element {
        public Simulation parentSimulation();
        public Class baseClass();
        public boolean wasAdded(); //indicates whether element was added to simulation (mediator)
        public void setAdded(boolean b);
        public void setName(String name);
        public String getName();
    }
    
    /**
     * Interface for a simulation element that can make a graphical component
     */
    public interface GraphicalElement extends Element {

        /**
         * Interface for a Simulation element that would be used in a simulation graphical user interface (GUI)
         * 
         * @param obj An object that might be used to specify the graphic that the GraphicalElement is to return.
         * In most cases the GraphicalElement ignores this parameter, and it can be set to null.
         * @return A Component that can be used in the GUI of a graphical simulation
         * @see Device
         * @see Display
         */
        public java.awt.Component graphic(Object obj);
    }//end of GraphicalElement

    public static final void makeAndDisplayFrame(Simulation sim) {
        javax.swing.JFrame f = new javax.swing.JFrame();
        f.setSize(600,500);
        f.getContentPane().add(sim.panel());
        f.pack();
        f.show();
        f.addWindowListener(Simulation.WINDOW_CLOSER);
    }
    
    public static final java.awt.event.WindowAdapter WINDOW_CLOSER 
        = new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        };
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        Default.ATOM_SIZE = 1.0;                   
	    IntegratorHard integratorHard = new IntegratorHard();
	    SpeciesDisks speciesDisks = new SpeciesDisks();
	    speciesDisks.setNMolecules(300);
	    Phase phase = new Phase();
	    Potential2 potential = new P2HardSphere();
	    Controller controller = new Controller();
	    DisplayPhase displayPhase = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard.new Timer(integratorHard.chronoMeter());
	    timer.setUpdateInterval(10);
        integratorHard.setTimeStep(0.01);
        displayPhase.setColorScheme(new ColorSchemeNull());
        for(Atom atom=phase.firstAtom(); atom!=null; atom=atom.nextAtom()) {
            atom.setColor(Constants.randomColor());
        }
        
        //this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		                                    
        Potential2.Agent potentialAgent = (Potential2.Agent)potential.getAgent(phase);
        potentialAgent.setIterator(new AtomPairIterator(phase));
		Simulation.instance.setBackground(java.awt.Color.yellow);
        Simulation.makeAndDisplayFrame(Simulation.instance);
        
     //   controller.start();
    }//end of main
}


