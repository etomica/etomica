package simulate;
import java.io.*;
import java.awt.*;
import java.beans.*;//for Graphics

public abstract class Species extends Component {

  String name;    // Name of species
  int speciesIndex;      // Numerical index of species
  protected int nMolecules;      //number of Molecules for this species
  public int nAtomsPerMolecule;
  int[] simulationPixelDimensions = {-1, -1}; // pixel width (0) and height (1) of simulation box less boundaries. 
  double[] simulationRunDimensions = {-1, -1}; //dimensions (0<...<=1.0) used in dynamics calculations.
  Phase parentPhase;
  Species nextSpecies, previousSpecies;
  Molecule[] molecule;
  Molecule firstMolecule, lastMolecule;  //first and last molecules
  Atom firstAtom, lastAtom;
  boolean fillVolume; // Fill volume when initializing if true
  double neighborUpdateSquareDisplacement = Double.MAX_VALUE;
  
  public Species() {
    this(20);  //default nMolecules
  }

  public Species(int n) {
    setDefaults();
    speciesIndex = 0;
    fillVolume = true;
    setNMolecules(n);    
  }

  //Define these methods in subclass to make and initialize appropriate type of molecule
  abstract void setDefaults();
  abstract void makeMolecules();
  abstract void initializeMolecules();
  
  public final int getNMolecules() {return nMolecules;}
  public final void setNMolecules(int n) {
    nMolecules = n;
    makeMolecules();
    orderMolecules();
    initializeMolecules();
  }
  
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
  }
        
        
  protected final void orderMolecules() {
    setFirstMolecule(molecule[0]);
    setLastMolecule(molecule[nMolecules-1]);
    for(int i=1; i<nMolecules; i++) {molecule[i-1].setNextMolecule(molecule[i]);}
  }
  
  private final void setFirstMolecule(Molecule m) {
    firstMolecule = m;
    firstAtom = (m != null) ? m.firstAtom : null;
  }
  
  private final void setLastMolecule(Molecule m) {
    lastMolecule = m;
    lastAtom = (m != null) ? m.lastAtom : null;
  }
  
  public final Species getNextSpecies() {return nextSpecies;}
  public final void setNextSpecies(Species s) {
    this.nextSpecies = s;
    s.previousSpecies = this;
    this.lastMolecule.setNextMolecule(s.firstMolecule);
  }
  public final Species getPreviousSpecies() {return previousSpecies;}
  
  public final Phase getPhase() {return parentPhase;}
    
  public final String getName() {return name;}
  public final void setName(String name) {this.name = name;}
  
  public final boolean isFillVolume() {return fillVolume;}
  public final void setFillVolume(boolean b) {fillVolume = b;}

  public final int getSpeciesIndex() {return speciesIndex;}
  public final void setSpeciesIndex(int index) {speciesIndex = index;}
  
  public final double getNeighborUpdateSquareDisplacement() {return neighborUpdateSquareDisplacement;}
  public final void setNeighborUpdateSquareDisplacement(double d) {neighborUpdateSquareDisplacement = d;}

  public double kineticEnergy() {
    double KE = 0;
    Molecule endMolecule = lastMolecule.getNextMolecule();
    for(Molecule m=firstMolecule; m!=endMolecule; m=m.getNextMolecule()) {
        KE += m.kineticEnergy();
    }
    return KE;
  }
  //paint should be executed only during design time, via call to paintComponents in Phase
  //At run time Phase calls draw instead
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
      
  // Design-time temporary workarounds
  protected double designTimeXDim = 1.0;  //duplicates element of space.dimensions[]
  protected double designTimeYDim = 1.0;

    public void setDesignTimeXDim(double x) {designTimeXDim = x;}
    public double getDesignTimeXDim() {return designTimeXDim;}
    public void setDesignTimeYDim(double y) {designTimeYDim = y;}
    public double getDesignTimeYDim() {return designTimeYDim;}
    
  
// Abstract Methods
  
  public abstract void draw(Graphics g, int[] origin, double scale);

  public abstract void initializeSpecies(Phase phase);
}