package simulate;
import java.awt.Container;
import java.awt.Graphics;

public class Phase extends Container {
        
    protected Space.Boundary boundary;
    private int iBoundary;
        
    public Phase() {
        iBoundary = 1;
        setLayout(null);
        setSize(300,300);
        atomCount = moleculeCount = 0;
        gravity = new Gravity(0.0);
        noGravity = true;
    }
    
    public void initialize(Simulation ps) {
        parentSimulation = ps;
        setBoundary(iBoundary);
        add(new ConfigurationSequential());  //default configuration
        potentialEnergy = new MeterPotentialEnergy();
        potentialEnergy.setUpdateInterval(Integer.MAX_VALUE);  //these meters are placed to permit phase to report its potential and kinetic energies
        kineticEnergy = new MeterKineticEnergy();
        kineticEnergy.setUpdateInterval(Integer.MAX_VALUE);    //change updateInterval if desired to use for averaging also
        add(potentialEnergy);
        add(kineticEnergy);
    }
        

    public final void setBoundary(int b) {
        iBoundary = b;
        if(parentSimulation != null) boundary = parentSimulation.space.makeBoundary(iBoundary);
    }
    public final int getBoundary() {return iBoundary;}
    public final Space.Boundary boundary() {return boundary;}
    
    public simulate.AtomPair makeAtomPair() {return makeAtomPair(null, null);}
    public simulate.AtomPair makeAtomPair(Atom a1, Atom a2) {return parentSimulation.space.makeAtomPair(boundary, a1, a2);}
    public final simulate.AtomPair.Iterator.A makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL) {return parentSimulation.space.makePairIteratorFull(boundary,iF,iL,oF,oL);}
    public final simulate.AtomPair.Iterator.A makePairIteratorHalf(Atom iL, Atom oF, Atom oL) {return parentSimulation.space.makePairIteratorHalf(boundary,iL,oF,oL);}
    public final simulate.AtomPair.Iterator.A makePairIteratorFull() {return parentSimulation.space.makePairIteratorFull(boundary);}
    public final simulate.AtomPair.Iterator.A makePairIteratorHalf() {return parentSimulation.space.makePairIteratorHalf(boundary);}
        
    public void paint(Graphics g, int[] origin, double scale) {} 
                    
    public final Space.Vector dimensions() {return boundary.dimensions();}
    public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
    public void inflate(double scale) {
        boundary.inflate(scale);
        for(Molecule m=firstMolecule(); m!=null; m=m.nextMolecule()) {
            m.inflate(scale);
        }
    }
    public void reflate(double scale) {
        boundary.inflate(1.0/scale);
        for(Molecule m=firstMolecule(); m!=null; m=m.nextMolecule()) {
            m.replace();
        }
    }
    
    public final Atom firstAtom() {
        Molecule m = firstMolecule();
        return (m != null) ? m.firstAtom() : null;
    }
    public final Atom lastAtom() {
        Molecule m = lastMolecule();
        return (m != null) ? m.lastAtom() : null;
    }
    public final Molecule firstMolecule() {
        for(Species.Agent s=firstSpecies; s!=null; s=s.nextSpecies()) {
            Molecule m = s.firstMolecule();
            if(m != null) {return m;}
        }
        return null;
    }
    public final Molecule lastMolecule() {
        for(Species.Agent s=lastSpecies; s!=null; s=s.previousSpecies()) {
            Molecule m = s.lastMolecule();
            if(m != null) {return m;}
        }
        return null;
    }
    public final Species.Agent firstSpecies() {return firstSpecies;}
    public final Species.Agent lastSpecies() {return lastSpecies;}
    public final Phase nextPhase() {return nextPhase;}
    public final Phase previousPhase() {return previousPhase;}
    /**
    * Sets the phase following this one in the linked list of phases.
    * Does not link species/molecule/atoms of the Phases
    *
    * @param p the phase to be designated as this phase nextPhase
    */
    public final void setNextPhase(Phase p) {
        this.nextPhase = p;
        p.previousPhase = this;
    }
          
    public final double getG() {return gravity.getG();}
    public void setG(double g) {
        gravity.setG(g);
        noGravity = (g == 0.0);
    }
          
    /**
    * Returns the temperature (in Kelvin) of this phase as computed via the equipartition
    * theorem from the kinetic energy summed over all (atomic) degrees of freedom
    */  
    public double kineticTemperature() {
        return (2./(double)(atomCount*parentSimulation.space.D()))*kineticEnergy.currentValue()*Constants.KE2T;
    }

    public void add(Configuration c){
        c.parentPhase = this;
        configuration = c;
        for(Species.Agent s=firstSpecies; s!=null; s=s.nextSpecies()) {
            configuration.add(s);
        }
    }
            
	public void add(Meter m) {
	    if(lastMeter != null) {lastMeter.setNextMeter(m);}
	    else {firstMeter = m;}
	    lastMeter = m;
	    meterCount++;
	    m.setPhase(this);
	    m.initialize();
	    if(parentSimulation != null && parentSimulation.haveIntegrator()) {
	        parentSimulation.controller.integrator.addIntegrationIntervalListener(m);
	    }
	}
        	
    public void add(Species.Agent species) {
//        species.configurationMolecule.initializeCoordinates();
        if(lastSpecies != null) {lastSpecies.setNextSpecies(species);}
        else {firstSpecies = species;}
        lastSpecies = species;
        for(Molecule m=species.firstMolecule(); m!=null; m=m.nextMolecule()) {moleculeCount++;}
        for(Atom a=species.firstAtom(); a!=null; a=a.nextAtom()) {atomCount++;}
        configuration.add(species);
    }
            
    public void add(Species s) {  //add species to phase if it doesn't appear in another phase
        s.parentSimulation = this.parentSimulation;
        Species.Agent agent = s.makeAgent(this);
        agent.setNMolecules(20);
        add(agent);
    }
                        
	// Returns ith meter in linked list of meters, with i=0 being the first meter
	public Meter getMeter(int i) {
	    if(i >= meterCount) {return null;}
	    Meter m = firstMeter;
        for(int j=i; --j>=0; ) {m = m.nextMeter();}  //get ith meter in list
        return m;
    }

    public void updateCounts() {
        moleculeCount = 0;
        atomCount = 0;
        for(Molecule m=firstMolecule(); m!=null; m=m.nextMolecule()) {moleculeCount++;}
        for(Atom a=firstAtom(); a!=null; a=a.nextAtom()) {atomCount++;}
    }
            
    /**
    * Object used to describe presence and magnitude of constant gravitational acceleration
    */
    public Gravity gravity;
    public boolean noGravity = true;
            
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
             
    Simulation parentSimulation;
    Meter firstMeter, lastMeter;
    private int meterCount = 0;
          
    public Configuration configuration;
          
    private Phase nextPhase;
    private Phase previousPhase;
          
    public Integrator integrator;
          
    public MeterPotentialEnergy potentialEnergy;
    public MeterKineticEnergy kineticEnergy;
 
}
        