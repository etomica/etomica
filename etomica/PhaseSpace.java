package simulate;
import java.awt.Container;
import java.awt.Graphics;

public abstract class PhaseSpace extends Container {

    public PhaseSpace() {
        setLayout(null);
        setSize(300,300);
        atomCount = moleculeCount = 0;
        gravity = new Gravity(0.0);
        noGravity = true;
        add(new ConfigurationSequential());  //default configuration
        potentialEnergy = new MeterPotentialEnergy();
        kineticEnergy = new MeterKineticEnergy();
        add(potentialEnergy);
    }
    
    public abstract int D();
    
    public abstract PhaseSpace.AtomCoordinate makeAtomCoordinate(Atom a);      //PhaseSpace prefix is redundant
    public abstract PhaseSpace.MoleculeCoordinate makeMoleculeCoordinate(Molecule m);
    public abstract PhaseSpace.Vector makeVector();
    public abstract simulate.AtomPair makeAtomPair(Atom a1, Atom a2);
    public abstract AtomPair.Iterator.A makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL);
    public abstract AtomPair.Iterator.A makePairIteratorHalf(Atom iL, Atom oF, Atom oL);
    public abstract AtomPair.Iterator.A makePairIteratorFull();
    public abstract AtomPair.Iterator.A makePairIteratorHalf();
    
    public abstract double volume();

//  Vector contains what is needed to describe a point in the space
    interface Vector {
        public void E(Vector u);
        public void E(double a);
        public void PE(Vector u);
        public void ME(Vector u);
        public void TE(Vector u);
        public void DE(Vector u);
        public void TE(double a);
        public void DE(double a);
        public double squared();
        public double dot(Vector u);
    }

//  Coordinate collects all vectors needed to describe point in phase space -- position and (maybe) momentum
    interface Coordinate {
        public void setNextCoordinate(Coordinate c);
        public void clearPreviousCoordinate();
        
        public void translateTo(Vector r);
        public void translateToward(Vector e, double distance);
        public void translateBy(Vector dr);
        public void displaceTo(Vector r);
        public void displaceBy(Vector dr);
        public void displaceWithin(double d);
        public void displaceToRandom();
        public void translateToRandom();
        public void randomizeMomentum(double temperature);
        public void replace();
        public void inflate(double s);
        public void accelerate(Vector dv);
        public double kineticEnergy();
        public Vector position();
        public Vector momentum();
        public Vector velocity();
    }
    
    interface AtomCoordinate extends Coordinate {      //cannot be a class here because must inherit from Coordinate as it is defined in the PhaseSpace subclass
        public Atom nextAtom();
        public Atom previousAtom();
        public Atom atom();
        public AtomCoordinate nextCoordinate();
        public AtomCoordinate previousCoordinate();
    }
    interface MoleculeCoordinate extends Coordinate {
        public Molecule nextMolecule();
        public Molecule previousMolecule();
        public Molecule molecule();
        public MoleculeCoordinate nextCoordinate();
        public MoleculeCoordinate previousCoordinate();
        public void displaceToRandom(Vector dim);
    }
        
 /**
  * Returns the temperature (in Kelvin) of this phase as computed via the equipartition
  * theorem from the kinetic energy summed over all (atomic) degrees of freedom
  */  
/*  public double getKineticTemperature() {   //move this functionality into a meter
    updateKineticEnergy();
    return (2./(double)(nAtomTotal*Space.D))*kineticEnergy*Constants.KE2T;
  }*/
  
  public final Atom firstAtom() {
     Molecule m = firstMolecule();
     return (m != null) ? m.firstAtom() : null;
  }
  public final Atom lastAtom() {
    Molecule m = lastMolecule();
    return (m != null) ? m.lastAtom() : null;
  }
  public final Molecule firstMolecule() {
    for(Species s=firstSpecies; s!=null; s=s.nextSpecies()) {
        Molecule m = s.firstMolecule();
        if(m != null) {return m;}
    }
    return null;
  }
  public final Molecule lastMolecule() {
    for(Species s=lastSpecies; s!=null; s=s.previousSpecies()) {
        Molecule m = s.lastMolecule();
        if(m != null) {return m;}
    }
    return null;
  }
  public final Species firstSpecies() {return firstSpecies;}
  public final Species lastSpecies() {return lastSpecies;}
  public final PhaseSpace nextPhaseSpace() {return nextPhaseSpace;}
  public final PhaseSpace previousPhaseSpace() {return previousPhaseSpace;}
 /**
  * Sets the phase following this one in the linked list of phases.
  * Does not link species/molecule/atoms of the Phases
  *
  * @param p the phase to be designated as this phase nextPhase
  */
  public final void setNextPhase(PhaseSpace p) {
    this.nextPhaseSpace = p;
    p.previousPhaseSpace = this;
  }
  
  public final double getG() {return gravity.getG();}
  public void setG(double g) {
    gravity.setG(g);
    noGravity = (g == 0.0);
  }
  
    public void add(Species species) {
        super.add(species);
        species.parentPhaseSpace = this;
        species.configurationMolecule.initializeCoordinates();
        configuration.add(species);
        if(lastSpecies != null) {lastSpecies.setNextSpecies(species);}
        else {firstSpecies = species;}
        lastSpecies = species;
        for(Molecule m=species.firstMolecule(); m!=null; m=m.nextMolecule()) {moleculeCount++;}
        for(Atom a=species.firstAtom(); a!=null; a=a.nextAtom()) {atomCount++;}
        if(species.getSpeciesIndex() > speciesCount-1) {setSpeciesCount(species.getSpeciesIndex()+1);}
    }
    
    /* Resizes potential arrays, keeping all elements already filled in, and
       setting to p1Null or p2IdealGas the newly formed elements
       */
    private void setSpeciesCount(int n) {
        Potential1 p1[] = new Potential1[n];
        Potential2 p2[][] = new Potential2[n][n];
        for(int i=0; i<speciesCount; i++) {
            p1[i] = potential1[i];
            p2[i][i] = potential2[i][i];
            for(int j=i+1; j<speciesCount; j++) {        //could use system arraycopy
                p2[i][j] = p2[j][i] = potential2[i][j];
                
            }
        }
        for(int i=speciesCount; i<n; i++) {
            p1[i] = p1Null;
            p2[i][i] = p2IdealGas;
            for(int j=0; j<n; j++) {
                p2[i][j] = p2[j][i] = p2IdealGas;
            }
        }
        potential1 = p1;
        potential2 = p2;
        speciesCount = n;
    }
  
    public void add(Potential1 p1) {
        super.add(p1);
        p1.setPhaseSpace(this);
        if(p1.speciesIndex+1 > speciesCount) {setSpeciesCount(p1.speciesIndex+1);}
        this.potential1[p1.speciesIndex] = p1;
    }
    
    public void add(Potential2 p2) {
        super.add(p2);
        p2.setPhaseSpace(this);
        int idx = Math.max(p2.species1Index,p2.species2Index);
        if(idx+1 > speciesCount) {setSpeciesCount(idx+1);}
        this.potential2[p2.species1Index][p2.species2Index] = p2;
        this.potential2[p2.species2Index][p2.species1Index] = p2;
    }
    
    public void add(Configuration c){
        c.parentPhaseSpace = this;
        configuration = c;
        for(Species s=firstSpecies; s!=null; s=s.nextSpecies()) {
            configuration.add(s);
        }
    }
    
	public void add(Meter m) {
	    if(lastMeter != null) {lastMeter.setNextMeter(m);}
	    else {firstMeter = m;}
	    lastMeter = m;
	    meterCount++;
	    m.phaseSpace = this;
	    m.initialize();
	    if(parentSimulation != null && parentSimulation.haveIntegrator()) {
	        parentSimulation.controller.integrator.addIntegrationIntervalListener(m);
	    }
	}
	
	// Returns ith meter in linked list of meters, with i=0 being the first meter
	public Meter getMeter(int i) {
	    if(i >= meterCount) {return null;}
	    Meter m = firstMeter;
        for(int j=i; --j>=0; ) {m = m.nextMeter();}  //get ith meter in list
        return m;
    }

    public void inflate(double scale) {
        dimensions.TE(scale);
    }

//    public final boolean isPeriodic() {return periodic;}  
        
   /**
    * @return the dimensions[] array
    */
    public Vector getDimensions() {return dimensions;}
   
   /** Set of vectors describing the displacements needed to translate the central image
    *  to all of the periodic images.  Returns a two dimensional array of doubles.  The
    *  first index specifies each perioidic image, while the second index indicates the
    *  x and y components of the translation vector.
    *  Likely to override in subclasses.
    *
    *  @param nShells the number of shells of images to be computed
    */
    public abstract double[][] imageOrigins(int nShells);
       
   /**
    * @return shift0
    * @see #shift0
    * @see Phase#paint
    */
    public double[][] getOverflowShifts(AtomCoordinate c, double distance) {return shift0;}  //called only if periodic
    
    public abstract void paint(Graphics g, int[] origin, double scale); 
    
 /**
  * Total number of species contained in this phase.
  */
  int speciesCount=0;
  
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
  * Object used to describe presence and magnitude of constant gravitational acceleration
  */
  public Gravity gravity;
  public boolean noGravity = true;
    
 /**
  * First species in the linked list of species in this phase.
  */
   private Species firstSpecies;
 
 /**
  * Last species in the linked list of species in this phase.
  */
  Species lastSpecies;
  
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
  
  private Potential1 p1Null = new P1Null();
  private Potential2 p2IdealGas = new P2IdealGas();
  public Configuration configuration;
  
  private PhaseSpace nextPhaseSpace;
  private PhaseSpace previousPhaseSpace;
  
  public Integrator integrator;
  
  public MeterPotentialEnergy potentialEnergy;
  public MeterKineticEnergy kineticEnergy;
    
 /**
  * Size of Phase (width, height) in Angstroms
  * Default value is 1.0 for each dimension.
  */
    protected Vector dimensions;

   /**
    * Volume of the phase, in Angstroms^D.
    *
    * @see #computeVolume
    */
    public double volume;
      
 /**
  * @see SpacePeriodicCubic#getOverflowShifts
  */
    protected final double[][] shift0 = new double[0][D()];
   
}