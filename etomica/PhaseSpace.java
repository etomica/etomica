package simulate;
import java.awt.Container;

public abstract class PhaseSpace extends Container {

    public PhaseSpace() {
        setLayout(null);
        setSize(300,300);
        setBackground(Color.white);
        atomCount = moleculeCount = 0;
        gravity = new Gravity(0.0);
        noGravity = true;
        add(new ConfigurationSequential());  //default configuration
        potentialEnergy = new MeterPotentialEnergy();
        add(potentialEnergy);
    }
    
    public abstract PhaseSpace.AtomCoordinate makeAtomCoordinate(Atom a);      //PhaseSpace prefix is redundant
    public abstract PhaseSpace.MoleculeCoordinate makeMoleculeCoordinate(Molecule m);
    public abstract PhaseSpace.AtomPair makeAtomPair(Atom a1, Atom a2);
    public abstract AtomPairIterator.A makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL);
    public abstract AtomPairIterator.A makePairIteratorHalf(Atom iL, Atom oF, Atom oL);
    public abstract AtomPairIterator.A makePairIteratorFull();
    public abstract AtomPairIterator.A makePairIteratorHalf();
    

//  Vector contains what is needed to describe a point in the space
    interface Vector {}

//  Coordinate collects all vectors needed to describe point in phase space -- position and (maybe) momentum
    interface Coordinate {
        public void setNextCoordinate(Coordinate c);
        public void clearPreviousCoordinate();
        
        public void translateTo(Vector r);
        public void translateBy(Vector dr);
        public void displaceTo(Vector r);
        public void displaceBy(Vector dr);
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
    }
    interface MoleculeCoordinate extends Coordinate {
        public Molecule nextMolecule();
        public Molecule previousMolecule();
        public void update();
    }
    
    public interface AtomPair {
        public double r2();
        public Atom atom1();
        public Atom atom2();
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
    for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
        Molecule m = s.firstMolecule();
        if(m != null) {return m;}
    }
    return null;
  }
  public final Molecule lastMolecule() {
    for(Species s=lastSpecies; s!=null; s=s.getPreviousSpecies()) {
        Molecule m = s.lastMolecule();
        if(m != null) {return m;}
    }
    return null;
  }
  public final Species firstSpecies() {return firstSpecies;}
  public final Species lastSpecies() {return lastSpecies;}
  public final Phase getNextPhase() {return nextPhase;}
  public final Phase getPreviousPhase() {return previousPhase;}
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
  
    public void add(Species species) {
        super.add(species);
        species.parentPhase = this;
        species.configurationMolecule.initializeCoordinates();
        configuration.add(species);
        if(lastSpecies != null) {lastSpecies.setNextSpecies(species);}
        else {firstSpecies = species;}
        lastSpecies = species;
        for(Molecule m=species.firstMolecule(); m!=null; m=m.getNextMolecule()) {moleculeCount++;}
        for(Atom a=species.firstAtom(); a!=null; a=a.getNextAtom()) {atomCount++;}
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
        p1.setPhase(this);
        if(p1.speciesIndex+1 > speciesCount) {setSpeciesCount(p1.speciesIndex+1);}
        this.potential1[p1.speciesIndex] = p1;
    }
    
    public void add(Potential2 p2) {
        super.add(p2);
        p2.setPhase(this);
        int idx = Math.max(p2.species1Index,p2.species2Index);
        if(idx+1 > speciesCount) {setSpeciesCount(idx+1);}
        this.potential2[p2.species1Index][p2.species2Index] = p2;
        this.potential2[p2.species2Index][p2.species1Index] = p2;
    }
    
    public void add(Configuration c){
        c.parentPhase = this;
        configuration = c;
        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
            configuration.add(s);
        }
    }
    
	public void add(Meter m) {
	    if(lastMeter != null) {lastMeter.setNextMeter(m);}
	    else {firstMeter = m;}
	    lastMeter = m;
	    meterCount++;
	    m.phase = this;
	    m.initialize();
	    if(parentSimulation != null && parentSimulation.haveIntegrator()) {
	        parentSimulation.controller.integrator.addIntegrationIntervalListener(m);
	    }
	}
	
	// Returns ith meter in linked list of meters, with i=0 being the first meter
	public Meter getMeter(int i) {
	    if(i >= meterCount) {return null;}
	    Meter m = firstMeter;
        for(int j=i; --j>=0; ) {m = m.getNextMeter();}  //get ith meter in list
        return m;
    }

    public final boolean isPeriodic() {return periodic;}  
        
   /**
    * @return the dimensions[] array
    */
    public Vector getDimensions() {return dimensions;}
   
   /**
    * Scales all dimensions by a constant multiplicative factor and recomputes volume
    *
    * @param scale the scaling factor. 
    */
    public void inflate(double scale) {
        Space.uTEa1(dimensions,scale);
        computeVolume();
    }
   
   /**
    * Computes the volume of the space based on the current values in
    * dimensions[], and assigns the result to <code>volume</code>.
    * Likely to override in subclasses.
    */
    protected void computeVolume() {
        volume = dimensions.x*dimensions.y;
    }
        
   /** Set of vectors describing the displacements needed to translate the central image
    *  to all of the periodic images.  Returns a two dimensional array of doubles.  The
    *  first index specifies each perioidic image, while the second index indicates the
    *  x and y components of the translation vector.
    *  Likely to override in subclasses.
    *
    *  @param nShells the number of shells of images to be computed
    */
    public double[][] imageOrigins(int nShells) {
        return new double[0][D];
    }
       
   /**
    * @return shift0
    * @see #shift0
    * @see Phase#paint
    */
    public Vector[][] getOverflowShifts(double[] r, double distance) {return shift0;}  //called only if periodic
    
   /**
    * Draws a light gray outline of the space if <code>visible</code> is set to
    * <code>true</code>.  The size of the outline is determined by 
    * <code>drawSize[]</code>.
    *
    * @param g      the graphics object to which the image is drawn
    * @param origin the coordinate origin (in pixels) for drawing the image
    * @see Phase#paint
    * @see #computeDrawSize
    */
    public void drawSpace(Graphics g, int[] origin, double scale) {
        g.setColor(Color.gray.brighter());
        g.drawRect(origin[0],origin[1],(int)(scale*dimensions.x)-1,(int)(scale*dimensions.y)-1);
    }
  
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
  
  private Phase nextPhase;
  private Phase previousPhase;
  
  public Integrator integrator;
  
  public MeterPotentialEnergy potentialEnergy;
    
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
      
}