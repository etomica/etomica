package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesCell extends Species {

  double[] speciesOrigin = {0,0}; //origin less depth of boundaries, shifts origin to accomodate boundaries.
  double mass;
  double diameter;
  double radius;
  Color color;
  
  public SpeciesDisk() {
    super();
    setDiameter(0.1);    
    setMass(1.0);
    setColor(Color.black);
  }
  
  public void setDefaults() {
    nTethersPerMolecule = 10;}
  
  //makeMolecules is called by constructor via setNMolecules
  void makeMolecules() {
    molecule = new Molecule[nMolecules];
    for(int i=0; i<nMolecules; i++) {molecule[i] = new MoleculeCell(this,nTethersPerMolecule);}
  }
  
  //initializeMolecules is called by constructor via setNMolecules
  void initializeMolecules() {
    initializeMolecules(diameter, mass, color);
  }
  void initializeMolecules(double d, double m, Color c) {
    setDiameter(d);    //call set methods to pass diameter and mass to atoms
    setMass(m);
    setColor(c);
  }
  
  // Exposed Properties
        
    public final Color getColor() {return color;}
    public final void setColor(Color c) {firstAtom.setColor(color);}
        
  public void draw(Graphics g, int[] origin, double scale) {
    double toPixels = scale*Phase.TO_PIXELS;
    int diameterP = (int)(toPixels*diameter);
    g.setColor(color);
    Atom nextSpeciesAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
        a.draw(g, origin, scale);
    }
  }

  //for default boundary widths of zero
  public void initializeCoordinates(double[] simulationRunDimensions) {
    double[] array = {0,0};
    initializeCoordinates(simulationRunDimensions, array);
  }

   // initializeSpecies called by this.paint at design time, and by phase.add at run time
  public void initializeSpecies(Phase phase) {
    parentPhase = phase;
    if(fillVolume) {
        setLocation(0,0);
        setSize(phase.getSize());
    }
    double[] simulationRunDimensions;
    if(Beans.isDesignTime()) {
        if(fillVolume) {
          speciesOrigin[0] = 0;
          speciesOrigin[1] = 0;
        }
        simulationRunDimensions = new double[Space.D];
        simulationRunDimensions[0] = designTimeXDim;
        simulationRunDimensions[1] = designTimeYDim;
    }
    else {
        speciesOrigin[0] = getLocation().x;
        speciesOrigin[1] = getLocation().y;
//        simulationRunDimensions = phase.space.getDimensions();
        simulationRunDimensions = new double[Space.D];
        simulationRunDimensions[0] = phase.space.getDimensions(0);
        simulationRunDimensions[1] = phase.space.getDimensions(1);
    }
    if(!fillVolume) {
        speciesOrigin[0] = getLocation().x/Phase.TO_PIXELS;
        speciesOrigin[1] = getLocation().y/Phase.TO_PIXELS;
        simulationRunDimensions[0] *= ((double)getSize().width)/(double)phase.getSize().width;
        simulationRunDimensions[1] *= ((double)getSize().height)/(double)phase.getSize().height;
    }
    initializeCoordinates(simulationRunDimensions);
  }

  public void initializeCoordinates(double[] simulationRunDimensions, double[] simulationOrigin) {
    double dx[] = new double[Space.D]
    dx[0] = simulationRunDimensions[0]/(double)nMolecules;
    dx[1] = 0.1;
    for(i=0; i<nMolecules; i++) {
        molecule[i].translate 
        molecule[i].atom[0].p[1] -= momentumSumY;
//        molecule[i].r[0] += 0.01*(2.0*rand.nextDouble()-1.0);
//        molecule[i].r[1] += 0.01*(2.0*rand.nextDouble()-1.0);
        }
  }
}


