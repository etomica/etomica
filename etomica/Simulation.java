//This class includes a main method to demonstrate its use
package etomica;

import javax.swing.*;
import java.util.*;
import etomica.units.UnitSystem;
import java.beans.Beans;

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
    
    public String getVersion() {return "Simulation:01.03.23";}
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
    private java.util.HashMap elementLists = new java.util.HashMap(16);
    
   /**
    * Object describing the nature of the physical space in which the simulation is performed
    */
    public Space space; //would like to make final, but compiler doesn't allow
    
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
     * List of all intra-molecular potentials that have been instantiated.
     */
    public LinkedList potential1List = new LinkedList();
    /**
     * List of all inter-molecular potentials that have been instantiated.
     */
    public LinkedList potential2List = new LinkedList();
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
	public final javax.swing.JPanel displayBoxPanel = new JPanel(new java.awt.GridLayout(0,1));
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
        elementLists.put(Potential1.class, potential1List);
        elementLists.put(Potential2.class, potential2List);
        elementLists.put(Species.class, speciesList);
        elementLists.put(Integrator.class, integratorList);
        elementLists.put(Phase.class, phaseList);
        elementLists.put(Controller.class, controllerList);
        elementLists.put(Display.class, displayList);
        elementLists.put(MeterAbstract.class, meterList);
        elementLists.put(Device.class, deviceList);
        p1Null = new P1Null(this);
        p2IdealGas = new P2IdealGas(this);
        elementCoordinator = new Mediator(this);
        setSize(800,550);
//        setLayout(new java.awt.FlowLayout());
        setLayout(new java.awt.BorderLayout());
        add(devicePanel, java.awt.BorderLayout.WEST);
        add(displayPanel, java.awt.BorderLayout.EAST);
        add(displayBoxPanel, java.awt.BorderLayout.CENTER);
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
     * @return the <code>nth</code> instantiated inter-molecular potential (indexing from zero)
     */
    public final Potential2 potential2(int n) {return (Potential2)potential2List.get(n);}
    /**
     * @return the <code>nth</code> instantiated intra-molecular potential (indexing from zero)
     */
    public final Potential1 potential1(int n) {return (Potential1)potential1List.get(n);}
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
    public final LinkedList potential1List() {return potential1List;}
    public final LinkedList potential2List() {return potential2List;}
    public final LinkedList displayList() {return displayList;}
    public final LinkedList deviceList() {return deviceList;}
  
    public void register(Simulation.Element element) {
        LinkedList list = (LinkedList)elementLists.get(element.baseClass());
        if(list.contains(element)) return;
        if(element instanceof P1Null) return;
        if(element instanceof P2IdealGas) return;
        element.setName(element.getClass().getName().substring(8) + Integer.toString(list.size()));
        list.add(element);
        allElements.add(element);
        
        if(element instanceof Species) {
            Species species = (Species)element;
            
            //update linked list of species
            if(lastSpecies != null) {lastSpecies.setNextSpecies(species);}
            else {firstSpecies = species;}
            lastSpecies = species;
            
            //set default index of species
            species.setSpeciesIndex(speciesCount()-1);
            
            //resize potential arrays
            int n = speciesCount()-1; //size of old arrayss
            Potential1[] new1 = new Potential1[n+1];
            Potential2[][] new2 = new Potential2[n+1][n+1];
            
            for(int i=0; i<n; i++) {
                new1[i] = potential1[i];
                new2[n][i] = p2IdealGas;
                new2[i][n] = p2IdealGas;
                for(int j=0; j<n; j++) {       
                    new2[i][j] = potential2[i][j];
                }
            }
            new1[n] = p1Null;
            new2[n][n] = p2IdealGas;
            potential1 = new1;
            potential2 = new2;
        }//end of if(element instanceof Species)
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
        
        //update linked list of species
        if(lastSpecies != null) {lastSpecies.setNextSpecies(species);}
        else {firstSpecies = species;}
        lastSpecies = species;
        
        //set default index of species
//        species.setSpeciesIndex(speciesCount()-1);
        
        /*if(Beans.isDesignTime()) {
            for(java.util.Iterator ip=Simulation.phaseList.iterator(); ip.hasNext(); ) {
                Phase p = (Phase)ip.next();
                p.remove(species.makeAgent(p));
            }
        }*/
    }
                
    /**
     * Method invoked in the constructor of a Potential1 object to list it with the simulation
     */
    public void unregister(Potential1 p1) {
        if(!potential1List.contains(p1)) return;
        if(p1 instanceof P1Null) return;
        int deleteIndex = potential1List.indexOf(p1);
        for (int i = deleteIndex; i < potential1List.size()-1; i++)
            potential1[i] = potential1[i+1];
        potential1List.remove(p1);
        allElements.remove(p1);
    }
                
    /**
     * Method invoked in the constructor of a Potential2 object to list it with the simulation
     */
    public void unregister(Potential2 p2) {
        if(!potential2List.contains(p2)) return;
        if(p2 instanceof P2IdealGas) return;
        potential2List.remove(p2);
        allElements.remove(p2);
    }
 
    public LinkedList allElements() {
        LinkedList list;
 //       synchronized(this) {
            list = (LinkedList)allElements.clone();
 //       }
        return list;
    }
            
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
    /**
     * Returns the potential governing the interaction of the given atoms
     */
    public final Potential getPotential(AtomPair pair) {
        Atom a1 = pair.atom1();
        Atom a2 = pair.atom2();
        return getPotential(a1,a2);
    }
    /**
     * Returns the potential governing the interaction of the given atoms
     */
    public final Potential getPotential(Atom a1, Atom a2) {
        if(a1 == a2) {
            return null;
//            return a1.parentPhase().potential();  //should rewrite AtomPair to hold phase
        }
        else if(a1.parentMolecule() == a2.parentMolecule()) {
            return potential1[a1.speciesIndex()].getPotential(a1,a2);}
        else {
            return potential2[a2.speciesIndex()][a1.speciesIndex()].getPotential(a1,a2);
        }
    }
    
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
    
    //Instances used as default molecular potentials, for molecules and pairs that don't have 
    //Potential1 and Potential2 set explicitly.
    private Potential1 p1Null;
    private P2IdealGas p2IdealGas;
    
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
    }

    /**
     * Method to construct a simple simulation of hard disks.
     * After this method is called, Simulation.instance will provide a handle to
     * the constructed simulation.  This feature is meant to show how a simulation
     * is constructed and accessed, and to provide a convenient instance to
     * aid testing during development
     *
     * @deprecated.  Use a preconstructed extension of Simulation instead.
     */
    public static void makeSimpleSimulation() {
	    IntegratorHard integratorHard1 = new IntegratorHard();
	    SpeciesDisks speciesDisks1 = new SpeciesDisks();
	    Phase phase1 = new Phase();
	    P2SimpleWrapper P2HardDisk1 = new P2SimpleWrapper(new PotentialHardDisk());
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
		Simulation.instance.setBackground(java.awt.Color.yellow);
    }//end of makeSimpleSimulation
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        Simulation.makeSimpleSimulation();  //for more general simulations, replace this call with
                                            //construction of the desired pieces of the simulation
                                            
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
}


