package etomica;

import java.awt.Graphics;
import java.awt.Color;
import java.util.Observable;
import java.util.Observer;
import java.util.LinkedList;
import java.util.Iterator;
import java.beans.Beans;
/**
 * A Phase collects all atoms that interact with one another; atoms in different phases do
 * not interact.  These are the important features of a Phase:<p> 
 * 
 * 1. It maintains an IteratorFactory that can be used to obtain iterators that loop over all
 * or some atoms and atom pairs in the phase.<p>
 * 
 * 2. It provides methods for addition and removal of molecules in the phase.<p>
 * 
 * 3. It holds a Boundary object, obtained from the governing Space, that defines the behavior
 * of atoms as they move into or across the boundary of the phase.<p>
 * 
 * 4. It holds a set of Meter objects that can be used to measure properties of that atoms
 * held by the phase.  An EnergyMeter is constructed by default to permit evaluation of
 * potential and kinetic energies of the atoms in the phase.<p>
 * 
 * 5. It has a Configuration object that determines the default initial configuration of the
 * atoms in the phase.<p>
 *
 * 6. It can contain one or more PotentialField objects that define one-body potentials affecting
 * the configurations of the atoms it contains.
 *
 * 
 * @author David Kofke
 * @see IteratorFactory
 * @see Space.Boundary
 * @see MeterAbstract
 * @see PotentialField
 */
public final class Phase implements Simulation.Element, Molecule.Container, java.io.Serializable {
        
    private Space.Boundary boundary;
    private transient final LinkedList meterList = new LinkedList();
    private PotentialField firstField = null;
    private PhaseAction.Inflate inflater;
    private String name;
    private final Simulation parentSimulation;
    private boolean added = false;
    
    public Phase() {
        this(Simulation.instance);
    }
    
    public Phase(Simulation sim) {
        parentSimulation = sim;
        inflater = new PhaseAction.Inflate(this);
        //don't use setIteratorFactory here to remove any possibility of complications with Observers
        if(sim.space() instanceof IteratorFactory.Maker) {
            iteratorFactory = ((IteratorFactory.Maker)sim.space()).makeIteratorFactory(this);
            if(iteratorFactory instanceof PotentialField.Maker) {
                addField(((PotentialField.Maker)iteratorFactory).makePotentialField(this));
            }
        }
        else {
            iteratorFactory = new IteratorFactory(this);
        }

        atomCount = moleculeCount = 0;

//        setBoundary(Space.Boundary.DEFAULT);
        setBoundary(parentSimulation().space().makeBoundary());

        if(parentSimulation.space().D() < 3) 
            add(new ConfigurationSequential());  //default configuration
        else
            add(new ConfigurationFcc());
        
        //Create and register energy meter; 
        //placed to permit phase to report its potential and kinetic energies
        energy = new MeterEnergy();
        energy.setActive(false);    //default has meter inactive, so it does not keep averages or respond to integrator interval events
        energy.setPhase(this);
        
        parentSimulation.register(this);
        
    }
    
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Phase.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
    public Space.Vector randomPosition() {return boundary.randomPosition();}
    
    public void setMoleculePosition(int i, Space.Coordinate r) {
        molecule(i).translateTo(r.position());
    }
    public void setMoleculePosition(Space.Coordinate[] r) {}
    /**
     * Returns a Space.Vector with the position of the <code>ith</code> molecule
     */
    public Space.Coordinate getMoleculePosition(int i) {
        return molecule(i).coordinate();}
    
    /**
     * Returns an array of Space.Vector with the center-of-mass position of each molecule
     * Not meant for computation-intensive use.
     */
    public Space.Coordinate[] getMoleculePosition() {
        Space.Coordinate[] positions = new Space.Coordinate[moleculeCount];
        Molecule molecule = firstMolecule();
        for(int i=0; i<moleculeCount; i++) {
            positions[i] = molecule.coordinate();
            molecule = molecule.nextMolecule();
        }
        return positions;
    }
    
    /**
     * Returns an array of the molecules currently in the phase, for those times
     * when you just need convenient indexed access to each molecule.
     * Use is discouraged.
     */
     public Molecule[] moleculeArray() {
        Molecule[] array = new Molecule[moleculeCount];
        int i=0;
        Molecule terminationMolecule = lastMolecule().nextMolecule();
        for(Molecule m=firstMolecule(); m!=terminationMolecule; m=m.nextMolecule()) {
            array[i++] = m;
        }
        return array;
     }
    
    /**
     * Returns the ith molecule in the linked list of molecules.
     * 0 returns the first molecule, and moleculeCount-1 returns the last.
     * An argument outside this range throws an IndexOutOfBoundsException
     */
    //could make more efficient by starting from first or last molecule, as appropriate
    public Molecule molecule(int i) {
        if(i >= moleculeCount || i < 0) 
            throw new IndexOutOfBoundsException("Index: "+i+
                                                ", Number of molecules: "+moleculeCount);
        Molecule m = firstMolecule();
        for(int count=i; count>0; count--) { 
            m = m.nextMolecule();
        }
        return m;
    }
    
    /**
     * Returns a randomly selected molecule from the phase.
     */
    public Molecule randomMolecule() {
        int i = (int)(moleculeCount * java.lang.Math.random());
        return molecule(i);
    }
        

    /**
     * Sets the boundary object of the phase.
     */
     public void setBoundary(Space.Boundary b) {
        boolean hasIntegrator = integrator != null;
        if(hasIntegrator) integrator.pause();
        if(boundary != null) { //remove any fields attributable to existing boundary
            for(PotentialField f=firstField; f!=null; f=f.nextField()) {
                if(f.maker()==boundary) {removeField(f);}
            }
        }
        boundaryMonitor.notifyObservers(b);
        boundary = b;
        boundary.setPhase(this);
        if(b instanceof PotentialField.Maker) {
            addField(((PotentialField.Maker)b).makePotentialField(this));
        }
        if(hasIntegrator) {
            integrator.reset();
            integrator.unPause();
        }
     }

    /**
     * Indicates the Boundary object of the phase
     * 
     * @return An integer that codes for the boundary via a static variable in Space.Boundary
     */
//    public final int getBoundary() {return iBoundary;}
    public final Space.Boundary getBoundary() {return boundary;}
    public final Space.Boundary boundary() {return boundary;}
    
    /**
     * Accessor of the phase's iterator factory object
     * 
     * @return The phase's iterator factory
     * @see IteratorFactory
     */
    public final IteratorFactory iteratorFactory() {return iteratorFactory;}

    /**
     * Sets the iterator factory of the phase
     * 
     * @param it
     * @see IteratorFactory
     */
    public final void setIteratorFactory(IteratorFactory it) {
        //notify observers before updating iteratorFactory in case observers want to close business with old iteratorFactory
        iteratorFactoryMonitor.notifyObservers(it);
        iteratorFactory = it;
        if(it instanceof PotentialField.Maker) {
            addField(((PotentialField.Maker)it).makePotentialField(this));
        }
    }
    
    /**
     * Accessor of the integrator governing the movement of the atoms in the phase.
     * 
     * @return The phase's integrator.
     * @see Integrator
     */
    public final Integrator integrator() {return integrator;}

    /**
     * Sets the integrator of the phase.
     * 
     * @param i The new integrator.
     */
    public void setIntegrator(Integrator i) {
           //need to define an exception for this
        if(!i.addPhase(this)) return; //addPhase will return false if the phase was not successfully added to the integrator
        //notify observers before updating integrator in case observers want to close business with old integrator
        integratorMonitor.notifyObservers(i);
//        integratorMonitor.fireEvent(new PhaseIntegratorEvent(i,PhaseIntegratorEvent.NEW_INTEGRATOR));
        if(integrator != null) integrator.removePhase(this);
        integrator = i;
        integrator.addIntervalListener(new PhaseAction.ImposePbc(this));
    }
                        
    public final Space.Vector dimensions() {return boundary.dimensions();}
    public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
    
    public void setDensity(double rho) {
        double vNew = moleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/parentSimulation().space().D());
        inflater.actionPerformed(scale);
    }
    public double getDensity() {return moleculeCount()/boundary.volume();}

    /**
     * @return the number of atoms in the Phase
     */
    public final int atomCount() {return atomCount;}
    
    /**
     * @return the number of molecules in the Phase
     */
    public final int moleculeCount() {return moleculeCount;}
    
    /**
     * @return the first atom in the linked list of atoms in this Phase
     */
    public final Atom firstAtom() {
        Molecule m = firstMolecule();
        return (m != null) ? m.firstAtom() : null;
    }
    
    /**
     * @return the last atom in the linked list of atoms in this Phase
     */
    public final Atom lastAtom() {
        Molecule m = lastMolecule();
        return (m != null) ? m.lastAtom() : null;
    }
    
    /**
     * @return the first molecule in the linked list of molecules in this Phase
     */
    public final Molecule firstMolecule() {
        for(Species.Agent s=firstSpecies; s!=null; s=s.nextSpecies()) {
            Molecule m = s.firstMolecule();
            if(m != null) {return m;}
        }
        return null;
    }
    
    /**
     * @return the last molecule in the linked list of molecules in this Phase
     */
    public final Molecule lastMolecule() {
        for(Species.Agent s=lastSpecies; s!=null; s=s.previousSpecies()) {
            Molecule m = s.lastMolecule();
            if(m != null) {return m;}
        }
        return null;
    }
    
    /**
     * @return the first Agent in the linked list of Species.Agents in this phase
     */
    public final Species.Agent firstSpecies() {return firstSpecies;}
    
    /**
     * @return the last Agent in the linked list of Species.Agents in this phase
     */
    public final Species.Agent lastSpecies() {return lastSpecies;}
          
    /**
    * Returns the temperature (in simulation units) of this phase as computed via the equipartition
    * theorem from the kinetic energy summed over all (atomic) degrees of freedom
    */  
    public double kineticTemperature() {
        return MeterTemperature.value(this);
    }
    
    public void setConfiguration(Configuration c) {
        configuration = c;
//        for(Species.Agent s=firstSpecies; s!=null; s=s.nextSpecies()) {
//            configuration.add(s);
//        }
        configuration.initializeCoordinates(this);
        iteratorFactory.reset();
    }
    
    public Configuration getConfiguration() {return configuration;}

    /**
     * @deprecated
     */
    public void add(Configuration c){setConfiguration(c);}
    
	/**
	 * Adds a meter to the list of meters working in this phase.
	 */
	public void addMeter(MeterAbstract m) {meterList.add(m);}

	/**
	 * Removes a meter from the list of meters working in this phase.
	 */
	public void removeMeter(MeterAbstract m) {meterList.remove(m);}

	/**
	 * Returns the list of meters working in this phase.
	 */
	public LinkedList getMeterList() {return meterList;}
    
    /**
     * Deploys the agent of a species in this phase
     */
    void addSpecies(Species.Agent species) {
        //set internal configuration of molecule
        if(species.parentSpecies().moleculeConfiguration != null) species.parentSpecies().moleculeConfiguration.initializeCoordinates(this);
        //add to linked list of species agents
        if(lastSpecies != null) {lastSpecies.setNextSpecies(species);}
        else {firstSpecies = species;}
        lastSpecies = species;
        //update molecule and atom counts
        updateCounts();
        //add species to configuration for this phase and notify iteratorFactory
        configuration.initializeCoordinates(this);
        iteratorFactory.reset();  
    }
    
    
    /**
     * Adds the given field to the phase.
     * If field is null, silently returns with no action taken.
     */
    public void addField(PotentialField field) {
        if(field == null) return;
        field.setNextField(firstField);
        firstField = field;
    }
    /**
     * Removes the given field from the phase.
     * If field is not in phase (or is null), silently returns with no action taken.
     * Value of nextField in removed field is kept intact, so field can be removed during a loop through
     * all fields without ruining loop.
     */
    public void removeField(PotentialField field) {
        if(firstField == null || field == null) return;
        if(field == firstField) {
            firstField = firstField.nextField();
        }
        else {
            PotentialField lastChecked = firstField;
            for(PotentialField f=firstField.nextField(); f!=null; f=f.nextField()) {
                if(f == field) { //remove f
                    lastChecked.setNextField(f.nextField());
                    break;
                }//end of if
                lastChecked = f;
            }//end of for
        }//end of else
    }//end of removeField
    
    /**
     * First field in a linked list of potential fields in this phase
     */
     public PotentialField firstField() {return firstField;}
    
    /**
     * Updates all the molecule and atoms counts for this phase.
     * Does not include wall atoms in the count.
     */
    public void updateCounts() {
        moleculeCount = 0;
        atomCount = 0;
        for(Molecule m=firstMolecule(); m!=null; m=m.nextMolecule()) {
            if(m.firstAtom.type instanceof AtomType.Wall) continue;
            moleculeCount++;
        }
        for(Atom a=firstAtom(); a!=null; a=a.nextAtom()) {
            if(a.type instanceof AtomType.Wall) continue;
            atomCount++;
        }
    }
    
    /**
     * Adds the given molecule to this phase, placing it in the molecule/atom linked lists
     * and removing it from the container it previously resided.
     * Looks up the molecule's species agent in this phase and passes it with the 
     * molecule to the two-argument addMolecule method.  If the agent is known already, 
     * that form of the method should be used instead.
     */
    public void addMolecule(Molecule m) {
        addMolecule(m, m.parentSpecies().getAgent(this));
    }
    
    /**
     * Adds the given molecule to this phase.
     * @param molecule the molecule to be added
     * @param s the species agent in this phase for the molecule's species.  If known, this agent can be provided to save the effort of looking it up.
     */
    public void addMolecule(Molecule m, Species.Agent s) {
        m.container().removeMolecule(m);
        m.setParentPhase(this);
        moleculeCount++;
        atomCount += m.atomCount;
        s.addMolecule(m);
        iteratorFactory.addMolecule(m);
    }
    
    /**
     * Removes molecule from phase.  
     * Should be called only by an addMolecule method of another container
     * Use deleteMolecule to remove molecule while not adding it to another phase (adds it to species reservoir)
     */
    public void removeMolecule(Molecule m) {
        removeMolecule(m, m.parentSpecies().getAgent(this));
    }
    public void removeMolecule(Molecule m, Species.Agent s) {
        m.setParentPhase(null);        
        moleculeCount--;
        atomCount -= m.atomCount;
        s.deleteMolecule(m);
        iteratorFactory.deleteMolecule(m);
    }
// deletes molecule by adding it to reservoir

public void deleteMolecule(Molecule m) {
        m.parentSpecies().reservoir().addMolecule(m);
    }
    /**
    * Synchronized version of deleteMolecule.  
    * Useful if molecules are being deleted by GUI events, rather than by integrator 
    */
    public final synchronized void deleteMoleculeSafely(Molecule m) {  //will this make deleteMolecule synchronized?
        deleteMolecule(m);
    }
    
    /**
    * Synchronized version of addMolecule
    * Useful if molecules are being added by GUI events, rather than by integrator 
    */
    public final synchronized void addMoleculeSafely(Molecule m) {
        addMolecule(m);
    }
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
    * First species in the linked list of species in this phase.
    */
    private Species.Agent firstSpecies;
         
    /**
    * Last species in the linked list of species in this phase.
    */
    Species.Agent lastSpecies;
          
    /**
    * Total number of atoms in this phase
    */
    public int atomCount;
         
    /**
    * Total number of molecules in this phase
    *
    * @see Species#addMolecule
    * @see Species#deleteMolecule
    */
    public int moleculeCount;
             
    public Configuration configuration;
          
    private Integrator integrator;
    
    private IteratorFactory iteratorFactory;
          
    public transient MeterEnergy energy;
        
    public Phase.Monitor integratorMonitor = new Phase.Monitor();
    public Phase.Monitor boundaryMonitor = new Phase.Monitor();
    public Phase.Monitor iteratorFactoryMonitor = new Phase.Monitor();
//    public SimulationEventManager integratorMonitor = new SimulationEventManager();
    
    public static class Monitor extends Observable implements java.io.Serializable {
        
        public void notifyObservers() {
            this.notifyObservers(null);
        }
        public void notifyObservers(Object obj) {
            setChanged();
            super.notifyObservers(obj);
        }
        public void addObserver(Observer o) {if(o != null) super.addObserver(o);}
    }
     
    /**
     * Class for constructing linked lists of Phases.
     * Each Linker points to one Phase and another Linker, the next one in the list.
     */
    public static class Linker implements java.io.Serializable {
        private final Phase phase;
        private Phase.Linker next = null;
        //Constructors
        public Linker(Phase a) {phase = a;}
        public Linker(Phase a, Linker l) {phase = a; next = l;}
        //Access methods
        public final Phase phase() {return phase;}
        public final Phase.Linker next() {return next;}
        public final void setNext(Phase.Linker l) {next = l;}
    }//end of Phase.Linker
    
    
} //end of Phase
        