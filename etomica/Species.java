package simulate;
import java.io.*;
import java.awt.*;
import java.beans.*;//for Graphics

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
    public void setNAtomsPerMolecule(int na)
    {
        nAtomsPerMolecule = na;
        setNMolecules(nMolecules);
    }
    
    public int getNAtomsPerMolecule() {return nAtomsPerMolecule;}
    
    public Species(int n, int na)
    {
        setSpeciesIndex(0);
        setFillVolume(true);
        nAtomsPerMolecule = na;
        setNMolecules(n);
        this.add(new ColorSchemeDefault(Color.red));
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
  protected int nAtomsPerMolecule;
  
 /**
  * Don't quite remember what this is used for; hope to deprecate it.
  */
  int[] simulationPixelDimensions = {-1, -1}; // pixel width (0) and height (1) of simulation box less boundaries. 
 
 /**
  * Don't quite remember what this is used for; hope to deprecate it.
  */
  double[] simulationRunDimensions = {-1, -1}; //dimensions (0<...<=1.0) used in dynamics calculations.
 
 /**
  * The phase in which this species resides.  Assigned in add method of Phase
  *
  */
  Phase parentPhase;
 
 /**
  * The Species following this one in the linked list of Species
  */
  Species nextSpecies;
 
 /**
  * The Species preceding this one in the linked list of Species
  */
  Species previousSpecies;
 
 /**
  * Array of Molecules collected by this Species.  This data structure
  * may be abandoned in future versions.
  */
  Molecule[] molecule;
 
 /**
  * First molecule in this Species
  *
  * @see Molecule
  */
  Molecule firstMolecule;
  
 /**
  * Last molecule in this Species
  *
  * @see Molecule
  */
  Molecule lastMolecule;
  
 /**
  * First Atom of the first Molecule in this Species
  *
  * @see Atom
  */
  public Atom firstAtom;
  
 /**
  * Last Atom of the last Molecule in this Species
  *
  * @see Atom
  */
  Atom lastAtom;
  
 /**
  * Use to determine initial placement of molecules.  If fillVolume is
  * true, initial coordinates fill all of Phase; otherwise initial-coordinate
  * cube may be shrunk at design time if desired.
  *
  * @see #draw
  */
  boolean fillVolume;
  
 /**
  * Parameter defining when neighbor list of a molecule of this species needs updating.
  * Current neighbor-list model is being deprecated.
  */
  double neighborUpdateSquareDisplacement = Double.MAX_VALUE;
  
  private transient final int[] shiftOrigin = new int[Space.D];     //work vector for drawing overflow images

 /**
  * Default constructor.  Creates species containing 20 molecules, each with 1 atom.
  */
  public Species() {
    this(20,1);
  }

 /**
  * Sets species defaults and creates molecules.  Called by default constructor.
  *
  * @param n number of molecules in species
  */
  public Species(int n) {
    this(n,1);
  }
  
  public void add(ColorScheme cs) {
    this.colorScheme = cs;
    for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {
        colorScheme.initializeAtomColor(a);
    }
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
  public final void setNMolecules(int n) {
    nMolecules = n;
    molecule = new Molecule[nMolecules];    //this array is to be deprecated
    for(int i=0; i<nMolecules; i++) {molecule[i] = makeMolecule();}
    orderMolecules();
    initializeMolecules();
  }
  
 /**
  * Creates new molecule from the Molecule class.  The number of atoms
  * per molecule is given by the current value of nAtomsPerMolecule.  
  * This class can be overridden to make molecules from a class other than
  * Molecule ("molecule" walls, for example).
  */
  public Molecule makeMolecule() {
    return new Molecule(this,nAtomsPerMolecule);
  }
 
 /**
  * Removes molecule from species, and updates atom and molecule linked lists.
  * Updates all values of first/last Molecule/Atom for species and
  * phase, if appropriate.  Also updates number-of-atom/molecule variables
  * for species and phase.
  * The Molecule array in species is not updated.
  * No measures are taken to remove this species if it holds zero molecules
  * after this molecule is deleted.
  *
  * @param m the molecule being deleted
  * @see #addMolecule
  */
  public final void deleteMolecule(Molecule m) {
    Molecule next = m.getNextMolecule();
    Molecule previous = m.getPreviousMolecule();
    if(m == firstMolecule) {
        if(nMolecules == 1) {setFirstMolecule(null);}
        else {setFirstMolecule(next);}
        if(m == parentPhase.firstMolecule) {
            parentPhase.setFirstMolecule();
            if(next != null) {
                next.previousMolecule = null;
                next.firstAtom.previousAtom = null;
            }
        }
    }
    if(m == lastMolecule) {
        if(nMolecules == 1) {setLastMolecule(null);}
        else {setLastMolecule(previous);}
        if(m == parentPhase.lastMolecule) {parentPhase.setLastMolecule();}
    }
    if(previous != null) {previous.setNextMolecule(next);}
    nMolecules--;
    parentPhase.nMoleculeTotal--;
    parentPhase.nAtomTotal -= m.nAtoms;
    m = null;
  }

 /**
  * Adds a molecule to this species and updates linked lists.  Does not handle
  * creation of molecule.  New molecule
  * becomes last molecule of species.  Updates first/last Molecule/Atom
  * for species and phase, if appropriate.
  * Not yet correctly implemented for use in a phase containing multiple
  * species (i.e., mixtures).  Adjusts total molecule and atom count in parent phase.
  *
  * @param m the molecule being added
  * @see deleteMolecule
  */
  public final void addMolecule(Molecule m) {
    if(nMolecules > 0) {
        m.setNextMolecule(lastMolecule.getNextMolecule());
        lastMolecule.setNextMolecule(m);
        setLastMolecule(m);
        if(parentPhase.lastMolecule == m.getPreviousMolecule()) {parentPhase.setLastMolecule();}
    }
    else {  //not suited for handling mixtures
        setFirstMolecule(m);
        setLastMolecule(m);
        parentPhase.setFirstMolecule();
        parentPhase.setLastMolecule();
    }
    nMolecules++;
    parentPhase.nMoleculeTotal++;
    parentPhase.nAtomTotal += m.nAtoms;
    colorScheme.initializeMoleculeColor(m);
  }
        
 /**
  * Link-list orders the molecules (and atoms) in this species.  Sets values
  * for first and last molecules.
  */ 
  protected final void orderMolecules() {
    setFirstMolecule(molecule[0]);
    setLastMolecule(molecule[nMolecules-1]);
    for(int i=1; i<nMolecules; i++) {molecule[i-1].setNextMolecule(molecule[i]);}
  }
  
 /**
  * Sets this species firstMolecule and firstAtom parameters
  *
  * @param m the molecule to be designated this species firstMolecule
  */
  private final void setFirstMolecule(Molecule m) {
    firstMolecule = m;
    firstAtom = (m != null) ? m.firstAtom : null;
  }
  
 /**
  * Sets this species lastMolecule and lastAtom parameters.
  *
  * @param m the molecule to be designated this species lastMolecule
  */
  private final void setLastMolecule(Molecule m) {
    lastMolecule = m;
    lastAtom = (m != null) ? m.lastAtom : null;
  }
  
 /**
  * @return the next species in the linked list of species.  Returns null if this is the last species.
  */
  public final Species getNextSpecies() {return nextSpecies;}
 
 /**
  * Sets the species following this one in the linked list of species.
  * Also links last Molecule/Atom of this species to the corresponding
  * first Molecule/Atom of the next species
  *
  * @param s the species to be designated as this species nextSpecies
  */
  public final void setNextSpecies(Species s) {
    this.nextSpecies = s;
    s.previousSpecies = this;
    this.lastMolecule.setNextMolecule(s.firstMolecule);
  }
 /**
  * @return the species preceding this one in the linked list of species.  Returns null if this is the first species.
  */
  public final Species getPreviousSpecies() {return previousSpecies;}
  
 /**
  * @return the phase in which this species resides
  */
  public final Phase getPhase() {return parentPhase;}
    
  public final String getName() {return name;}
  public final void setName(String name) {this.name = name;}
  
  public final boolean isFillVolume() {return fillVolume;}
  public final void setFillVolume(boolean b) {fillVolume = b;}

  public final int getSpeciesIndex() {return speciesIndex;}
  public final void setSpeciesIndex(int index) {speciesIndex = index;}
  
  public final double getNeighborUpdateSquareDisplacement() {return neighborUpdateSquareDisplacement;}
  public final void setNeighborUpdateSquareDisplacement(double d) {neighborUpdateSquareDisplacement = d;}

 /**
  * Computes and returns the kinetic energy summed over all molecules in this
  * species
  *
  * @return  kinetic energy in (amu)(Angstrom)<sup>2</sup>(ps)<sup>-2</sup>
  * @see Molecule#kineticEnergy
  */
  public double kineticEnergy() {
    double KE = 0;
    Molecule endMolecule = lastMolecule.getNextMolecule();
    for(Molecule m=firstMolecule; m!=endMolecule; m=m.getNextMolecule()) {
        KE += m.kineticEnergy();
    }
    return KE;
  }
  
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
        if(getParent() != null) {
            Phase par = (Phase)getParent();
            initializeSpecies(par);
            int[] origin = new int[Space.D];
            double scale = par.getNominalScale()/(2.0*par.getImageShells()+1);
            origin[0] = (par.getSize().width - Phase.toPixels(scale*designTimeXDim))/2;
            origin[1] = (par.getSize().height - Phase.toPixels(scale*designTimeYDim))/2;
            draw(g,origin,scale);
        }
    }
  }
      
  /**
   * Design-time variable
   *
   * @see #paint
   */
  protected double designTimeXDim = 1.0;
  
  /**
   * Design-time variable
   *
   * @see #paint
   */
  protected double designTimeYDim = 1.0;

    public void setDesignTimeXDim(double x) {designTimeXDim = x;}
    public double getDesignTimeXDim() {return designTimeXDim;}
    public void setDesignTimeYDim(double y) {designTimeYDim = y;}
    public double getDesignTimeYDim() {return designTimeYDim;}
  
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
    Atom nextSpeciesAtom = lastAtom.getNextAtom();
    Molecule last = lastMolecule.getNextMolecule();
    for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
        colorScheme.setAtomColor(a);
        a.draw(g,origin,scale);
        /*
        int xP = origin[0] + (int)(toPixels*(a.r[0]-radius));
        int yP = origin[1] + (int)(toPixels*(a.r[1]-radius));
        g.fillOval(xP,yP,diameterP,diameterP);
        */
    }
    if(parentPhase.drawOverflowImages) {
        for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
            double[][] shifts = parentPhase.space.getOverflowShifts(a.r,a.radius);
            for(int i=0; i<shifts.length; i++) {
                /*
               int xP = origin[0] + (int)(toPixels*(shifts[i][0]+a.r[0]-radius));
               int yP = origin[1] + (int)(toPixels*(shifts[i][1]+a.r[1]-radius));
               g.fillOval(xP,yP,diameterP,diameterP);
               */
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
  public abstract void initializeSpecies(Phase phase);
  
  public void initializeSpecies() {
    if(parentPhase != null) initializeSpecies(parentPhase);
  }
  
  public void setBounds(int x, int y, int width, int height) {
    Rectangle r = getBounds();
    if(r.x!=x || r.y!=y || r.width!=width || r.height!=height) {  
        super.setBounds(x, y, width, height);
        initializeSpecies();
    }
  }
}