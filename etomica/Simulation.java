//This class includes a main method to demonstrate its use
package simulate;

import javax.swing.*;
import java.util.*;
import simulate.units.UnitSystem;
import java.beans.Beans;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

/**
 * The main class that organizes the elements of a molecular simulation.
 * Contains many static fields and methods that provide a common point of reference
 * for all simulation objects.  Holds a single space object that is referenced in
 * many places to obtain spatial elements such as vectors and boundaries.  Also
 * holds an object that specifies the unit system used to default all I/O.  A single
 * instance of Simulation is held as a static field.  This instance is used to hold
 * graphic objects that might be displayed in a GUI.  These may be shown by adding 
 * this instance to a Container object.
 *
 * @author David Kofke
 */
public class Simulation extends javax.swing.JPanel implements java.io.Serializable {
    
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
    
   /**
    * Object describing the nature of the physical space in which the simulation is performed
    */
    public Space space; //would like to make final, but compiler doesn't allow
    
    /**
     * List of all controllers that have been instantiated.
     */
    public static LinkedList controllerList = new LinkedList();
    /**
     * List of all phases that have been instantiated.
     */
    public static LinkedList phaseList = new LinkedList();
    /**
     * List of all species that have been instantiated.
     */
    public static LinkedList speciesList = new LinkedList();
    /**
     * List of all displays that have been instantiated.
     */
    public static LinkedList displayList = new LinkedList();
    /**
     * List of all devices that have been instantiated.
     */
    public static LinkedList deviceList = new LinkedList();
    /**
     * List of all integrators that have been instantiated.
     */
    public static LinkedList integratorList = new LinkedList();
    /**
     * List of all intra-molecular potentials that have been instantiated.
     */
    public static LinkedList potential1List = new LinkedList();
    /**
     * List of all inter-molecular potentials that have been instantiated.
     */
    public static LinkedList potential2List = new LinkedList();
    /**
     * List of all meters that have been instantiated.
     */
    public static LinkedList meterList = new LinkedList();
    /**
     * List of all graphicalElements (Devices or Displays) that have been instantiated.
     */
    public static LinkedList graphicalElementList = new LinkedList();
    /**
     * List of all simulation elements.
     */
     private static LinkedList allElements = new LinkedList();

    //default unit system for I/O (internal calculations are all done in simulation units)
    private static UnitSystem unitSystem = new UnitSystem.Sim();
    
	public final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
	public final javax.swing.JPanel displayBoxPanel = new JPanel(new java.awt.GridLayout(0,1));
    public final javax.swing.JPanel devicePanel = new JPanel(new java.awt.GridLayout(0,1));
    

    /**
     * A static instance of a Simulation, suitable as a default Container of graphical simulation elements.
     * The add method of Simulation (which is inherited from Panel) is not
     * static, but it can be invoked for this static instance of Simulation.  In this way
     * graphical elements can be collected in this instance by default.  This instance
     * can then be added to an applet or application.
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
        name = "Simulation" + Integer.toString(instanceCount);
        instanceCount++;
        p1Null = new P1Null(this);
        p2IdealGas = new P2IdealGas(this);
//        elementCoordinator = new MediatorOneIntegrator(this);
        elementCoordinator = new Mediator(this);
        setSize(800,550);
//        setLayout(new GridLayout(1,2));
        setLayout(new java.awt.FlowLayout());
        add(devicePanel);
        add(displayPanel);
        add(displayBoxPanel);
        //workaround for JTabbedPane bug in JDK 1.2
        displayPanel.addChangeListener(
            new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent event) {
                    displayPanel.validate();
                }
        });
    }
    
    /**
     * Sets the type of space used in <b>all</b> simulations
     */
//    public void setSpace(Space s) {
//        space = s;
//    }
    
    private void writeObject(ObjectOutputStream out) throws java.io.IOException {
        out.defaultWriteObject();
//        out.writeObject(space);  //removed when space was made non-static
        out.writeObject(phaseList);
        out.writeObject(integratorList);
        out.writeObject(speciesList);
        out.writeObject(potential1List);
        out.writeObject(potential2List);
        out.writeObject(controllerList);
        out.writeObject(displayList);
        out.writeObject(meterList);
        out.writeObject(deviceList);
        out.writeObject(graphicalElementList);
 //       out.writeObject(potential1);
 //       out.writeObject(potential2);
        out.writeObject(firstSpecies);
        out.writeObject(lastSpecies);
 //       out.writeObject(p1Null);
 //       out.writeObject(p2IdealGas);
    }
    
    private void readObject(ObjectInputStream in) throws java.io.IOException, java.lang.ClassNotFoundException {
        in.defaultReadObject();
//        space = (Space)in.readObject();
        phaseList = (LinkedList)in.readObject();
        integratorList = (LinkedList)in.readObject();
        speciesList = (LinkedList)in.readObject();
        potential1List = (LinkedList)in.readObject();
        potential2List = (LinkedList)in.readObject();
        controllerList = (LinkedList)in.readObject();
        displayList = (LinkedList)in.readObject();
        meterList = (LinkedList)in.readObject();
        deviceList = (LinkedList)in.readObject();
        graphicalElementList = (LinkedList)in.readObject();
 //       potential1 = (Potential1[])in.readObject();
 //       potential2 = (Potential2[][])in.readObject();
        firstSpecies = (Species)in.readObject();
        lastSpecies = (Species)in.readObject();
 //       p1Null = (Potential1)in.readObject();
 //       p2IdealGas = (P2SimpleWrapper)in.readObject();
    }
    
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
    public static final Phase phase(int n) {return (Phase)phaseList.get(n);}
    /**
     * @return the <code>nth</code> instantiated species (indexing from zero)
     */
    public static final Species species(int n) {return (Species)speciesList.get(n);}
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
    public static final Controller controller(int n) {return (Controller)controllerList.get(n);}
    /**
     * @return the <code>nth</code> instantiated integrator (indexing from zero)
     */
    public static final Integrator integrator(int n) {return (Integrator)integratorList.get(n);}
    /**
     * @return the <code>nth</code> instantiated meter (indexing from zero)
     */
    public static final MeterAbstract meter(int n) {return (MeterAbstract)meterList.get(n);}
    /**
     * @return the <code>nth</code> instantiated display (indexing from zero)
     */
    public static final Display display(int n) {return (Display)displayList.get(n);}
    /**
     * @return the <code>nth</code> instantiated device (indexing from zero)
     */
    public static final Device device(int n) {return (Device)deviceList.get(n);}
  
    /**
     * Method invoked in the constructor of a Controller object to list it with the simulation
     */
    public static void register(Controller c) {
        if(controllerList.contains(c)) return;
        controllerList.add(c);
        allElements.add(c);
    }
                
    /**
     * Method invoked in the constructor of a Display object to list it with the simulation
     */
    public static void register(Display d) {
        if(displayList.contains(d)) return;
        displayList.add(d);
        allElements.add(d);
        graphicalElementList.add(d);
    }
    
    /**
     * Method invoked in the constructor of a Device object to list it with the simulation
     */
    public static void register(Device d) {
        if(deviceList.contains(d)) return;
        deviceList.add(d);
        allElements.add(d);
        graphicalElementList.add(d);
    }
              
    /**
     * Method invoked in the constructor of a Phase object to list it with the simulation
     */
    public static void register(Phase p) {
        if(phaseList.contains(p)) return;
        phaseList.add(p);
        allElements.add(p);
        
/*        if(Beans.isDesignTime()) {
            for(java.util.Iterator is=Simulation.speciesList.iterator(); is.hasNext(); ) {
                Species species = (Species)is.next();
                p.addSpecies(species.makeAgent(p));
            }
        }*/
    }
    
    /**
     * Method invoked in the constructor of an Integrator object to list it with the simulation
     */
    public static void register(Integrator i) {
        if(integratorList.contains(i)) return;
        integratorList.add(i); 
        allElements.add(i);
    }
    
    /**
     * Method invoked in the constructor of a MeterAbstract object to list it with the simulation
     */
    public static void register(MeterAbstract m) {
        if(meterList.contains(m)) return;
        meterList.add(m); 
        allElements.add(m);
    }
                  
    /**
     * Method invoked in the constructor of a Species object to list it with the simulation
     */
    public void register(Species species) {
        if(speciesList.contains(species)) return;
        //register with manager
        speciesList.add(species);
        allElements.add(species);
        
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
        
/*        if(Beans.isDesignTime()) {
            for(java.util.Iterator ip=Simulation.phaseList.iterator(); ip.hasNext(); ) {
                Phase p = (Phase)ip.next();
                p.addSpecies(species.makeAgent(p));
            }
        }*/
        potential1 = new1;
        potential2 = new2;
    }
                
    /**
     * Method invoked in the constructor of a Potential1 object to list it with the simulation
     */
    public static void register(Potential1 p1) {
        if(potential1List.contains(p1)) return;
        if(p1 instanceof P1Null) return;
        potential1List.add(p1);
        allElements.add(p1);
    }
                
    /**
     * Method invoked in the constructor of a Potential2 object to list it with the simulation
     */
    public static void register(Potential2 p2) {
        if(potential2List.contains(p2)) return;
        if(p2 instanceof P2IdealGas) return;  //to get an ideal-gas potential registered, make one using P2SimpleWrapper
        potential2List.add(p2);
        allElements.add(p2);
    }
                
/*****
 *
 * Beginning of Bryan's Remove methods
 *
 */    
 
    /**
     * Method invoked in the constructor of a Controller object to list it with the simulation
     */
    public void unregister(Controller c) {
        if(!controllerList.contains(c)) return;
        controllerList.remove(c);
        allElements.remove(c);
        
        /*for (int i = 0; i < Simulation.instance.getComponentCount(); i++) {
            System.out.println(Simulation.instance.getComponent(i).getName());
	        if (Simulation.instance.getComponent(i).getName() == "displayPanel"){
	            System.out.println("eat one");
                javax.swing.JPanel dp = ((javax.swing.JPanel)Simulation.instance.getComponent(i));
                dp.remove(dp.getComponent(0));
                dp.repaint();
                elementCoordinator.completed = false;
            }
        }*/
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
        graphicalElementList.remove(d);
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
        graphicalElementList.remove(d);
    }
              
    /**
     * Method invoked in the constructor of a Phase object to list it with the simulation
     */
    public void unregister(Phase p) {
        if(!phaseList.contains(p)) return;
        phaseList.remove(p);
        allElements.remove(p);
        
        /*if(Beans.isDesignTime()) {
            for(java.util.Iterator is=Simulation.speciesList.iterator(); is.hasNext(); ) {
                Species species = (Species)is.next();
                p.add(species.makeAgent(p));
            }
        }*/
    }
    
    /**
     * Method invoked in the constructor of an Integrator object to list it with the simulation
     */
    public void unregister(Integrator i) {
        if(!integratorList.contains(i)) return;
        integratorList.remove(i); 
        allElements.remove(i);
    }
    
    /**
     * Method invoked in the constructor of a MeterAbstract object to list it with the simulation
     */
    public void unregister(MeterAbstract m) {
        if(!meterList.contains(m)) return;
        meterList.remove(m); 
        allElements.remove(m);
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
 
/*****
 *
 * Endof Bryan's Remove methods
 *
 */                
    public static LinkedList allElements() {
        LinkedList list;
 //       synchronized(this) {
            list = (LinkedList)allElements.clone();
 //       }
        return list;
    }
            
    /**
     * Accessor of the first species in the linked list of species
     */
    public static final Species firstSpecies() {return firstSpecies;}
    /**
     * Accessor of the last species in the linked list of species
     */
    public static final Species lastSpecies() {return lastSpecies;}
    /**
     * @return the number of species that have been constructed and registered with the simulation
     */
     public static final int speciesCount() {return speciesList.size();}
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
    public static Species firstSpecies;
     
    /**
    * Last species in the linked list of species in this phase.
    */
    public static Species lastSpecies;
    
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
    * A marker interface to indicate that the class is an element of a Simulation
    */
    public interface Element {
        public Simulation parentSimulation();
        public Class baseClass();
        public boolean wasAdded(); //indicates whether element was added to simulation (mediator)
        public void setAdded(boolean b);
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


