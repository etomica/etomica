package simulate;
import java.io.*;
import java.awt.*;
import java.beans.*;//for Graphics
import java.util.Random;

/**
 * Each Phase contains one or more Species.
 * A Species is a collection of identically formed Molecules.  Subclasses of
 * Species differ in the type of molecules that they collect.  Subclassing
 * Species is the primary way to define new types of molecules (Molecule is
 * subclassed only to define molecules formed from unusual (e.g., anisotropic)
 * atoms).
 * The number of Molecules in a Species may be changed at run time.
 * One or more Species are collected together to form a Phase.
 * Interactions among all molecules in a Phase are defined by associating an
 * intermolecular potential to pairs of Species (or to a single Species to
 * define interactions among its molecules). Intramolecular potentials are
 * also associated with a Species, and are used to describe intra-molecular 
 * atomic interactions.
 *
 * Species are beans.
 *
 * @author David Kofke
 * @author C. Daniel Barnes
 *
 * @see Molecule
 * @see Potential1
 * @see Potential2
 * @see Phase
 */
 
public abstract class Species extends Container {


//No base-class constructor
  
  public void add(ColorScheme cs) {
    this.colorScheme = cs;
    for(Atom a=firstAtom(); a!=terminationAtom(); a=a.nextAtom()) {
        colorScheme.initializeAtomColor(a);
    }
  }
  
  public void add(ConfigurationMolecule cm) {
    cm.parentSpecies = this;
    this.configurationMolecule = cm;
    cm.initializeCoordinates();
  }
  
 /**
  * Copies all values for atoms sizes, masses, colors, etc., to the corresponding
  * molecule and species parameters.
  * This method is called by setNMolecules once the molecules have been made 
  * and linked-list ordered.
  */
  abstract void initializeMolecules();
 
 /**
  * @return the number of molecules for this species
  */
  public final int getNMolecules() {return nMolecules;}
  
 /**
  * Sets the number of molecules for this species.  Makes the given number
  * of new molecules, linked-list orders and initializes them.
  * Any previously existing molecules for this species are abandoned, as are
  * their links to molecules of next or previous species.
  *
  * @param n  the new number of molecules for this species
  * @see #makeMolecule
  * @see #orderMolecules
  * @see #initializeMolecules
  * @see #deleteMolecule
  * @see #addMolecule
  */
  public void setNMolecules(int n) {
    nMolecules = n;
    if(nMolecules == 0) {
        firstMolecule = null;
        lastMolecule = null;
        return;
    }
    firstMolecule = makeMolecule();
    lastMolecule = firstMolecule;
    for(int i=1; i<nMolecules; i++) {
        lastMolecule.setNextMolecule(makeMolecule());
        lastMolecule = lastMolecule.nextMolecule();
    }
    spareMolecule = makeMolecule();
    initializeMolecules();
  }
  
 /**
  * Creates new molecule of this species.  
  * If spareMolecule is not null, it is returned; else a new molecule is created
  * The coordinates of the molecule cannot be assumed to have any particular form
  * (e.g., the molecule is not necessarily at the origin)
  */
  public abstract Molecule makeMolecule();
/*    if(spareMolecule == null) {
        return new Molecule(this, atomsPerMolecule);
    }
    else {
        Molecule m = spareMolecule;
        spareMolecule = null;
        return m;
    }
 */ 
  
  /**
   * Chooses a molecule randomly from Species
   *
   * @return the randomly seleted molecule
   */
  public final Molecule randomMolecule() {
     int i = (int)(rand.nextDouble()*nMolecules);
     Molecule m = firstMolecule;
     for(int j=i; --j>=0; ) {m = m.nextMolecule();}
     return m;
  }
    
  /**
   * Chooses a molecule randomly from Species, and deletes it (removes it from the linked list)
   *
   * @return the deleted molecule
   */
  public final Molecule deleteMolecule() {
    Molecule m = this.randomMolecule();
    deleteMolecule(m);
    return m;
  }
    
  /**
   * Synchronized version of deleteMolecule.  
   * Useful if molecules are being deleted by GUI events, rather than by integrator 
   */
  public final synchronized void deleteMoleculeSafely(Molecule m) {  //will this make deleteMolecule synchronized?
    deleteMolecule(m);
  }
  
 /**
  * Removes molecule from species, and updates atom and molecule linked lists.
  * Updates all values of first/last Molecule/Atom for species and
  * phase, if appropriate.  Also updates number-of-atom/molecule variables
  * for species and phase.
  * No measures are taken to remove this species if it holds zero molecules
  * after this molecule is deleted.
  *
  * @param m the molecule being deleted
  * @see #addMolecule
  */
  
  public final void deleteMolecule(Molecule m) {
    if(m.parentSpecies != this) {
        System.out.println("Error:  attempt to delete molecule from incorrect species");
        return;
    }
    Molecule next = m.nextMolecule();
    Molecule previous = m.previousMolecule();
    if(m == firstMolecule) {
        if(nMolecules == 1) {firstMolecule = null;}  //deleting the first and only molecule of the species
        else {firstMolecule = next;}                 //deleting first molecule, but others are present
    }
    if(m == lastMolecule) {
        if(nMolecules == 1) {lastMolecule = null;}
        else {lastMolecule = previous;}
    }
    if(previous != null) {previous.setNextMolecule(next);} //reconnect linked list if not at beginning
    else if(next != null) {next.clearPreviousMolecule();}  //beginning of list; no previous molecule for next
    nMolecules--;
    m.parentSpecies = null;        //line deleted because of spareMolecule
    m.setNextMolecule(null);
    m.clearPreviousMolecule();
    parentPhaseSpace.moleculeCount--;
    parentPhaseSpace.atomCount -= m.nAtoms;
//    if(spareMolecule == null) spareMolecule = m;
  }

  /**
   * Synchronized version of addMolecule
   * Useful if molecules are being added by GUI events, rather than by integrator 
   */
  public final synchronized void addMoleculeSafely(Molecule m) {
    addMolecule(m);
  }
  
 /**
  * Adds a molecule to this species and updates linked lists.  Does not handle
  * creation of molecule.  New molecule
  * becomes last molecule of species.  Updates first/last Molecule
  * for species, if appropriate.
  * Not yet correctly implemented for use in a phase containing multiple
  * species (i.e., mixtures).  Adjusts total molecule and atom count in parent phase.
  *
  * @param m the molecule being added
  * @see deleteMolecule
  */
  public final void addMolecule(Molecule m) {
    if(nMolecules > 0) {
        m.setNextMolecule(lastMolecule.nextMolecule());
        lastMolecule.setNextMolecule(m);
        lastMolecule = m;
    }
    else {  //m is going to be the only molecule in species
        firstMolecule = m;
        lastMolecule = m;
        m.setNextMolecule(null); 
        for(Species s=this.nextSpecies(); s!=null; s=s.nextSpecies()) { //loop forward in species, looking for next molecule
            if(s.firstMolecule() != null) {
                m.setNextMolecule(s.firstMolecule());
                break;
            }
        }
        m.clearPreviousMolecule();
        for(Species s=this.previousSpecies(); s!=null; s=s.previousSpecies()) { //loop backward in species, looking for previous molecule
            if(s.lastMolecule() != null) {
                s.lastMolecule.setNextMolecule(m);
                break;
            }
        }
    }
    nMolecules++;
    m.parentSpecies = this;
    parentPhaseSpace.moleculeCount++;
    parentPhaseSpace.atomCount += m.nAtoms;
    initializeMolecules();
    colorScheme.initializeMoleculeColor(m);
  }
  
 /**
  * Makes a new molecule, assigns it a random position and velocity, and adds it to the species.  
  *
  * @return the added molecule
  */
  public Molecule addMolecule() {
    Molecule m = makeMolecule();
    configurationMolecule.initializeCoordinates(m);   //initialize internal coordinates
    m.setCOM(parentPhaseSpace.randomPosition());       //place at random position
    parentPhaseSpace.configuration.initializeMomentum(m);  //initialize momentum
    addMolecule(m);
    return m;
  }
        
 /**
  * @return the next species in the linked list of species.  Returns null if this is the last species.
  */
  public final Species nextSpecies() {return nextSpecies;}
 
 /**
  * Sets the species following this one in the linked list of species.
  * Also links last Molecule/Atom of this species to the corresponding
  * first Molecule/Atom of the next species
  *
  * @param s the species to be designated as this species nextSpecies
  */
  public final void setNextSpecies(Species s) {
    this.nextSpecies = s;
    Molecule last = lastMolecule();
    if(s==null) {
        if(last!=null) last.setNextMolecule(null); 
        return;
    }
    s.previousSpecies = this;
    if(last != null) {last.setNextMolecule(s.firstMolecule);}
  }
 /**
  * @return the species preceding this one in the linked list of species.  Returns null if this is the first species.
  */
  public final Species previousSpecies() {return previousSpecies;}
  
 /**
  * @return the phase in which this species resides
  */
  public final PhaseSpace parentPhaseSpace() {return parentPhaseSpace;}
    
  public final String getName() {return name;}
  public final void setName(String name) {this.name = name;}
  
//  public final boolean isFillVolume() {return fillVolume;}
//  public final void setFillVolume(boolean b) {fillVolume = b;}

  public final int getSpeciesIndex() {return speciesIndex;}
  public final void setSpeciesIndex(int index) {speciesIndex = index;}
  
 /**
  * Method used to draw molecules of species to screen at design time.
  * Invoked via call to paintComponents in Phase
  * At run time Phase calls draw instead
  * Uses designTimeXDim and designTimeYDim to determine values that at runtime
  * are given by Space.dimensions[0 or 1].  Also takes the value of Space.scale
  * to be unity.  This is done because at design
  * time Space is not accessible to Species.  One hopes that newer versions of
  * the Java language (or improvements in our understanding of it) 
  * will permit a more elegant means to handle this matter.
  *
  * @param g Graphics object to which molecule images are drawn
  * @see #draw
  */
  public void paint(Graphics g) {
    if(Beans.isDesignTime()) {
        if(parentPhaseSpace != null) {
            int[] origin = new int[Simulation.D];
            double scale = 1.0;
            origin[0] = (parentPhaseSpace.getSize().width - Phase.toPixels(scale*1.0))/2;  //1.0 -> dimension.x
            origin[1] = (parentPhaseSpace.getSize().height - Phase.toPixels(scale*1.0))/2; //       dimension.y
            draw(g,origin,scale);
        }
    }
  }
   
  /**
   * Draws all molecules of the species using current values of their positions.
   *
   * @param g         graphics object to which molecules are drawn
   * @param origin    origin of drawing coordinates (pixels)
   * @param scale     factor determining size of drawn image relative to
   *                  nominal drawing size
   * @see Atom#draw
   * @see Molecule#draw
   * @see Phase#paint
   */
  
  public void draw(Graphics g, int[] origin, double scale) {

    double toPixels = scale*Phase.TO_PIXELS;
 /*   
    int diameterP = (int)(toPixels*diameter);
    g.setColor(color);
    */
    Atom nextSpeciesAtom = lastAtom().nextAtom();
    Molecule last = lastMolecule.nextMolecule();
    for(Atom a=firstAtom(); a!=nextSpeciesAtom; a=a.nextAtom()) {
        colorScheme.setAtomColor(a);
        a.draw(g,origin,scale);
        /*
        int xP = origin[0] + (int)(toPixels*(a.r[0]-radius));
        int yP = origin[1] + (int)(toPixels*(a.r[1]-radius));
        g.fillOval(xP,yP,diameterP,diameterP);
        */
    }
    if(DisplayConfiguration.DRAW_OVERFLOW && (atomGenerator instanceof AtomDisk)) {
        for(Atom a=firstAtom(); a!=nextSpeciesAtom; a=a.nextAtom()) {
            double[][] shifts = parentPhaseSpace.getOverflowShifts(a.r,((AtomDisk)a).getRadius());  //should instead of radius have a size for all AtomC types
            for(int i=0; i<shifts.length; i++) {
               shiftOrigin[0] = origin[0] + (int)(toPixels*shifts[i][0]);
               shiftOrigin[1] = origin[1] + (int)(toPixels*shifts[i][1]);
               a.draw(g,shiftOrigin,scale);
            }
        }
    } 
  }
  
 /**
  * Method that sets the initial coordinates of the molecules in the
  * species.  If fillVolume is <code>true</code>, will set coordinate
  * origin and value to fill entire space in Phase; if fillVolume is <code>
  * false</code>, will not do this, and thereby permit adjustment of 
  * origin and size of occupied volume at design time.
  * Called by paint at design time, and by Phase.add at run time, and by this.setBounds (override of superclass)
  */
//  public abstract void initializeSpecies(Phase phase);
  
  public void setBounds(int x, int y, int width, int height) {
///    Rectangle r = getBounds();
///    if(r.x!=x || r.y!=y || r.width!=width || r.height!=height) {  
///        super.setBounds(x, y, width, height);
//        if(parentPhase != null) initializeSpecies(parentPhase);
///    }
  }
 /* 
  public void addNotify() {
    super.addNotify();
    if(getParent() instanceof Phase) {
       parentPhase = (Phase)getParent();
    }
    else {  //use an exception here?
        System.out.println("Error:  Species can be added only to a Phase");
    }
  }
*/  
    public void setNAtomsPerMolecule(int na)
    {
        atomsPerMolecule = na;
        setNMolecules(nMolecules);
    }
    
    public int getAtomsPerMolecule() {return atomsPerMolecule;}
    
  
  public final Molecule firstMolecule() {return firstMolecule;}
  public final Molecule lastMolecule() {return lastMolecule;}
  
  /**
   * Used to terminate loops over molecules in species
   */
  public final Molecule terminationMolecule() {
    return (lastMolecule == null) ? null : lastMolecule.nextMolecule();
  }  
  public final Atom firstAtom() { //return firstAtom;
    return (firstMolecule == null) ? null : firstMolecule.firstAtom();
  }
  public final Atom lastAtom() { //return lastAtom;
    return (lastMolecule == null) ? null : lastMolecule.lastAtom();
  }
  
  /**
   * Used to terminate loops over atoms in species
   */
  public final Atom terminationAtom() {
    Atom last = lastAtom();
    return (last == null) ? null : last.nextAtom();
  }

   /**
  * A name to be associated with the species.  Use is optional.
  */
  String name;
  
  public ColorScheme colorScheme;
 
 /**
  * A unique integer that identifies the species.  Used to associate
  * inter- and intra-molecular potentials to the Species.
  * Default value is often 0, but may use setDefaults method in subclass
  * to override this to another value if convenient; may also be set 
  * at design time
  * No two Species should have a common speciesIndex, but no check is made
  * to prevent this. Failure to properly assign values will cause mis-association
  * of Potentials and Species
  * 
  * @see Phase#potential2
  * @see Phase#potential1
  */
  int speciesIndex;
  
 /**
  * Number of molecules for this species
  */
  protected int nMolecules;
  
 /**
  * Number of atoms in each molecule of this species
  */
  protected int atomsPerMolecule;
      
 /**
  * The phase in which this species resides.  Assigned in add method of Phase
  *
  */
  PhaseSpace parentPhaseSpace;
 
 /**
  * The Species following this one in the linked list of Species
  */
  Species nextSpecies;
 
 /**
  * The Species preceding this one in the linked list of Species
  */
  Species previousSpecies;
  
 /**
  * First molecule in this Species
  *
  * @see Molecule
  */
  protected Molecule firstMolecule;
  
 /**
  * Last molecule in this Species
  *
  * @see Molecule
  */
  protected Molecule lastMolecule;
  
 /**
  * Extra molecule for doing test insertions and facilitating molecule exchanges
  */
  public Molecule spareMolecule;   //no longer used
      
  private transient final int[] shiftOrigin = new int[Simulation.D];     //work vector for drawing overflow images

 /**
  * Object responsible for setting default configuration of atoms in molecule
  */
  public ConfigurationMolecule configurationMolecule;
  
  private final Random rand = new Random();  //for choosing molecule at random

}