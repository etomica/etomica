package simulate;
import java.awt.*;//Because this class subclasses Component.
import java.io.*;//for Serializable and Stream.
import java.util.*; //for type Vector
import java.beans.*;//for PropertyChangeSupport and other beans convenience classes.
import java.util.*; //for neighborList
import java.awt.event.*;

/**
 * A Phase comprises a Space, one or more Species, and Potentials that characterize
 * how the species interact.  It may be viewed as a data structure that organizes
 * this information and which may be acted upon by an Integrator (which is the object
 * that contains information about how to do dynamics and ensemble sampling).
 * The Phase is responsible for knowing how to report instantaneous values (as opposed
 * to averages) of various mechanical and structural quantities (such as its energy
 * or a histogram of atom distances), and it knows how to draw itself to the screen.
 *
 * A Phase is a bean.  It is placed inside a Simulation, and it contains one Space, 
 * one or more Species, and one or more Potential1 and Potential2 objects.
 *
 * @author David Kofke
 * @author C. Daniel Barnes
 *
 * @see Simulation
 * @see Space
 * @see Species
 * @see Potential1
 * @see Potential2
 */
public final class Phase extends Container {

 /**
  * Conversion factor from simulation unit of length (Angstroms) to screen pixels (default = 300).
  * It is a class variable, so any instance of Phase that changes its value will
  * affect how all Phases are drawn to the screen.  The instance variable maySetTO_PIXELS
  * determines whether a phase has permission to change TO_PIXELS.
  */
  public static double TO_PIXELS = 300.;
 
 /**
  * Flag specifying whether phase computes energies, etc., using neighbor lists.  
  * Neighbor list feature is under revision, so presently their use should be avoided.
  */
  public  boolean useNeighborList;
 
 /**
  * Flag indicating whether this phase has permission to change the value of TO_PIXELS
  * Default value is <code>true</code>.
  *
  * @see #TO_PIXELS
  */
  private boolean maySetTO_PIXELS = true;
  
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
   
  public boolean updatedForces, updatedPotentialEnergy, updatedKineticEnergy;
  public boolean updatedNeighbors, updatedFutureNeighbors;
  
//  private double potentialEnergy
  private double kineticEnergy;  // must access with get methods
 
 /**
  * The Space object associated with this Phase.
  */
  Space space;
 
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
  public int nAtomTotal;
 
 /**
  * Total number of molecules in this phase
  *
  * @see Species#addMolecule
  * @see Species#deleteMolecule
  */
  public int nMoleculeTotal;
     
  Simulation parentSimulation;
  Meter firstMeter, lastMeter;
  private int nMeters = 0;
  
  private Potential1 p1Null = new P1Null();
  private Potential2 p2IdealGas = new P2IdealGas();
  public Configuration configuration;
  
  private Phase nextPhase;
  private Phase previousPhase;
  
  public Integrator integrator;
  
  public MeterPotentialEnergy potentialEnergy;
    
  public Phase() {
    setLayout(null);
    setSize(300,300);
    setBackground(Color.white);
    updatedForces = updatedPotentialEnergy = updatedKineticEnergy = false;
    updatedNeighbors = updatedFutureNeighbors = false;
    useNeighborList = false;
    nAtomTotal = nMoleculeTotal = 0;
    gravity = new Gravity(0.0);
    noGravity = true;
    add(new ConfigurationSequential());  //default configuration
    System.out.println("adding potentialenergy");
    potentialEnergy = new MeterPotentialEnergy();
    add(potentialEnergy);
  }
  
 /**
  * Returns the temperature (in Kelvin) of this phase as computed via the equipartition
  * theorem from the kinetic energy summed over all (atomic) degrees of freedom
  */  
  public double getKineticTemperature() {
    updateKineticEnergy();
    return (2./(double)(nAtomTotal*Space.D))*kineticEnergy*Constants.KE2T;
  }
  
  public final SpaceAtom firstAtom() {
     Molecule m = firstMolecule();
     return (m != null) ? m.firstAtom() : null;
  }
  public final SpaceAtom lastAtom() {
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
//    this.lastSpecies().setNextSpecies(p.firstSpecies);
  }
  
  /** 
   * Resets the class variable Phase.TO_PIXELS so that the x-dimension of this Phase's
   * space exactly fills this Phase's drawing region when space is unscaled
   * setTO_PIXELS = (width of this Phase, in pixels)/(width of Space, in Angstroms)
   *
   * @see #paint
   */
  public void resetTO_PIXELS() {
    if(maySetTO_PIXELS && space!=null) {
        setTO_PIXELS((double)getSize().width/space.getDimensions(0));
    }
  }
  
  /**
   * Sets the value of maySetTO_PIXELS for this phase.  If set to true,
   * resetTO_PIXELS is then called
   *
   * @param b the new value of maySetTO_PIXELS
   * @see #TO_PIXELS
   * @see #resetTO_PIXELS
   */
  public void setMaySetTO_PIXELS(boolean b) {
    maySetTO_PIXELS = b;
    resetTO_PIXELS();
  }
  
 /**
  * @return current value of maySetTO_PIXELS
  * @see #TO_PIXELS
  */
  public boolean getMaySetTO_PIXELS() {return maySetTO_PIXELS;}
 
 /**
  * Sets TO_PIXELS to the given value if maySetTO_PIXELS is <code>true</code>
  */
  public void setTO_PIXELS(double t) {if(maySetTO_PIXELS) {TO_PIXELS = t;}}
 
 /**
  * @return the current value of TO_PIXELS
  */
  public double getTO_PIXELS() {return TO_PIXELS;}
 
 /**
  * Performs a conversion from Angstroms to pixels
  * @param x a length, in simulation units (Angstroms)
  * @return x*TO_PIXELS, the input, converted to pixels
  */
  public static int toPixels(double x) {return (int)(TO_PIXELS*x);}

  public boolean getUseNeighborList() {return useNeighborList;}
  public void setUseNeighborList(boolean b) {useNeighborList = b;}
  
  public final double getG() {return gravity.getG();}
  public void setG(double g) {
    gravity.setG(g);
    noGravity = (g == 0.0);
  }
  
  // STOPPED commenting here
  
  
  // end of size-change method overrides

    public void add(Species species) {
        super.add(species);
        species.parentPhase = this;
        species.configurationMolecule.initializeCoordinates();
//        if(space != null) species.initializeSpecies(this);   delete if works ok
        configuration.add(species);
        if(lastSpecies != null) {lastSpecies.setNextSpecies(species);}
        else {firstSpecies = species;}
        lastSpecies = species;
        for(Molecule m=species.firstMolecule(); m!=null; m=m.getNextMolecule()) {nMoleculeTotal++;}
        for(Atom a=species.firstAtom(); a!=null; a=a.getNextAtom()) {nAtomTotal++;}
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
            for(int j=i+1; j<speciesCount; j++) {
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
  
     public void add(Space space) {
        super.add(space);
        this.space = space;
        space.setParentPhase(this);
        int w = getSize().width;
        int h = getSize().height;
        resetTO_PIXELS();
//        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {s.initializeSpecies(this);}
        for(int i = 0; i < speciesCount; i++) {
            if(potential1[i] != null) potential1[i].setSpace(space);
            for(int j = i; j < speciesCount; j++) {
                if(potential2[i][j] != null) potential2[i][j].setSpace(space);
            }
        }   
    }
    
    public void add(Potential1 p1) {
        super.add(p1);
        p1.setPhase(this);
        if(space != null) p1.setSpace(space);
        if(p1.speciesIndex+1 > speciesCount) {setSpeciesCount(p1.speciesIndex+1);}
        this.potential1[p1.speciesIndex] = p1;
    }
    
    public void add(Potential2 p2) {
        super.add(p2);
        p2.setPhase(this);
        if(space != null) p2.setSpace(space);
        int idx = Math.max(p2.species1Index,p2.species2Index);
        if(idx+1 > speciesCount) {setSpeciesCount(idx+1);}
        this.potential2[p2.species1Index][p2.species2Index] = p2;
        this.potential2[p2.species2Index][p2.species1Index] = p2;
/*        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
            if(p2.species1Index == s.getSpeciesIndex() || p2.species2Index == s.getSpeciesIndex()) {
                double d = 0.25*p2.skinThickness*p2.skinThickness;
                if(d < s.getNeighborUpdateSquareDisplacement()) {
                    s.setNeighborUpdateSquareDisplacement(d);
                }
            }
          }
        */
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
	    nMeters++;
	    m.phase = this;
	    m.initialize();
	    if(parentSimulation != null && parentSimulation.haveIntegrator()) {
	        parentSimulation.controller.integrator.addIntegrationIntervalListener(m);
	    }
	}
	
	// Returns ith meter in linked list of meters, with i=0 being the first meter
	public Meter getMeter(int i) {
	    if(i >= nMeters) {return null;}
	    Meter m = firstMeter;
        for(int j=i; --j>=0; ) {m = m.getNextMeter();}  //get ith meter in list
        return m;
    }

    // Updates neighbor list for all molecules
    // Only preceding molecules are included in a given molecule's neighbor list
    public void updateNeighbors() {
        for(Molecule m1=firstMolecule(); m1!=null; m1=m1.getNextMolecule()) {
            m1.clearNeighborList();
            for(Molecule m2=firstMolecule(); m2!=m1; m2=m2.getNextMolecule()) {
                if(potential2[m1.getSpeciesIndex()][m2.getSpeciesIndex()].isNeighbor(m1,m2)) {
                    m1.addNeighbor(m2);
                }
            }
        }
        updatedNeighbors = true;
    }
        
    public void updateForces() {
/*      if(updatedForces) {return;}
      potentialEnergy = 0.0;
      for(Molecule m1=firstElement; m1!=null; m1=m1.getNextMolecule()) {
        m1.zeroForce();
        Molecule m2;
//        for(Molecule m2=firstElement; m2!=m1; m2=m2.getNext()) {
        for(Enumeration enum=m1.getNeighborList(); enum.hasMoreElements();) {
           m2 = (Molecule)enum.nextElement();
           PairInteraction pair = potential[m1.getSpeciesIndex()][m2.getSpeciesIndex()].computePairInteraction(m1,m2);
           m2.addForce(pair.force);
           m1.subtractForce(pair.force);
           potentialEnergy += pair.energy;
        }
      }
      updatedForces = true;
      updatedPotentialEnergy = true;
*/    }
    
    public void updateKineticEnergy() {
      kineticEnergy = 0.0;
      for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
        kineticEnergy += s.kineticEnergy();
      }
      updatedKineticEnergy = true;
    }
 
 //  This method to be eliminated
 // Works only for one species
/*    public boolean overlap(Atom atom, double energy) {
      energy = 0.0;
      for(Atom a=firstAtom(); a!=atom; a=a.getNextAtom()) {
         double u = 0.0;
         Potential pot = potential2[0][0].getPotential(atom,a);
         if(pot.overlap(atom,a,u)) {return true;}
         else {energy += u;}
      }
      for(Atom a=atom.getNextAtom(); a!=null; a=a.getNextAtom()) {
         double u = 0.0;
         Potential pot = potential2[0][0].getPotential(atom,a);
         if(pot.overlap(atom,a)) {return true;}
         else {energy += u;}
      }      
      return false;
    }
*/        
    public void updatePotentialEnergy() {
/*      potentialEnergy = 0.0;
      for(Molecule m1=firstMolecule; m1!=null; m1=m1.getNextMolecule()) {
        for(Molecule m2=firstMolecule; m2!=m1; m2=m2.getNextMolecule()) {
          PairInteraction pair = potential[m1.getSpeciesIndex()][m2.getSpeciesIndex()].computePairInteraction(m1,m2);
          potentialEnergy += pair.energy;
        }
      }
      updatedPotentialEnergy = true;
*/    }

//needs work for efficiency and multiatomics   
/*    public double potentialEnergy() {  
        double energy = 0.0;
        for(Molecule m1=firstMolecule(); m1!=null; m1=m1.getNextMolecule()) {
           energy += m1.potentialEnergy();
        }
        return (0.5*energy);
      }
        */

    public double getKineticEnergy() {
        if(!updatedKineticEnergy) {updateKineticEnergy();}
        return kineticEnergy;
    }
 /*   
    public double getPotentialEnergy() {
        if(!updatedPotentialEnergy) {updatePotentialEnergy();}
        return potentialEnergy;
    }
    
    public double getTotalEnergy() {
        return getKineticEnergy() + getPotentialEnergy();
    }
      */
}
    