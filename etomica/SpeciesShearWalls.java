package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesShearWalls extends Species {

  double velocity = 50.0;
  double[] speciesOrigin = {0.,0.}; //origin less depth of boundaries, shifts origin to accomodate boundaries.
  double[] mass;      //make arrays when doublearrayeditor in place
  double[] diameter;
  double[] radius;
  Color color;
  
  double L;    //atom-atom separation
  public int nAtomsPerMolecule;
 
  public SpeciesShearWalls() {
    super(2);
    setColor(Color.gray);
    L = 0.20;               //default bondlength
    name = "Shear walls";
  }
  
  public void setDefaults() {
    nAtomsPerMolecule = 5;
    mass = new double[nAtomsPerMolecule];
    diameter = new double[nAtomsPerMolecule];
    radius = new double[nAtomsPerMolecule];
    for(int i=0; i<nAtomsPerMolecule; i++) {
        mass[i] = 1.0e20;
        diameter[i] = 0.1;
    }
  }
  
  //makeMolecules is called by constructor via setNMolecules
  void makeMolecules() {
    molecule = new Molecule[nMolecules];
    for(int i=0; i<nMolecules; i++) {molecule[i] = new Molecule(this,nAtomsPerMolecule);}
  }
  
  //initializeMolecules is called by constructor via setNMolecules
  void initializeMolecules() {
    initializeMolecules(diameter, mass, color);
  }
  void initializeMolecules(double[] d, double[] m, Color c) {
    setDiameter(d);    //call set methods to pass diameter and mass to atoms
    setMass(m);
    setColor(c);
  }
    
  // Exposed Properties
  
  public final double getVelocity() {return velocity;}
  public final void setVelocity(double v) {velocity = v;}
  
/*  public final double getMass() {return mass;}
  public final void setMass(double mass) {
    this.mass = mass;
    for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {a.setMass(mass);}
  }
  */  
  public final double[] getDiameter() {return diameter;}
  public final double getDiameter(int i) {return diameter[i];}
  public final void setDiameter(double[] d) {for(int i=0; i<nAtomsPerMolecule; i++) {setDiameter(i,d[i]);}}
  public final void setDiameter(int i, double value) {
    diameter[i] = value;
    radius[i] = 0.5*value;
    for(Molecule m=firstMolecule; m!=lastMolecule.getNextMolecule(); m=m.getNextMolecule()) {
        m.atom[i].setDiameter(value);
    }
  }
        
  public final Color getColor() {return color;}
  public final void setColor(Color c) {
    color = c;
    for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {a.setColor(c);}
  }

  public double[] getSpeciesOrigin() {return speciesOrigin;}
  public void setSpeciesOrigin(double[] speciesOrigin) {
    this.speciesOrigin = speciesOrigin;
  }
  public double getSpeciesOrigin(int i) {return speciesOrigin[i];}
  public void setSpeciesOrigin(int i, double value) {
    this.speciesOrigin[i] = value;
  }
        
  public double[] getMass() {return mass;}
  public double getMass(int i) {return mass[i];}
  public void setMass(double[] mass) {for(int i=0; i<nAtomsPerMolecule; i++) setMass(i,mass[i]);}
  public void setMass(int i, double value) {
    mass[i] = value;
    for(Molecule m=firstMolecule; m!=lastMolecule.getNextMolecule(); m=m.getNextMolecule()) {
        m.atom[i].setMass(mass[i]);
        m.updateMass();
    }        
  }
 
  public double getBondLength() {return L;}
  public void setBondLength(double L) {this.L = L;}

  public final double kineticEnergy() {return 0.0;}
  
  public void draw(Graphics g, int[] origin, double scale) {
    double toPixels = scale*Phase.TO_PIXELS;
    g.setColor(color);
    int o0 = origin[0];  //copy to scalars to avoid array-bounds check in loop
    int o1 = origin[1];
    Atom endAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=endAtom; a=a.getNextAtom()) {
        int xP = o0 + (int)(toPixels*(a.r[0]-a.radius));
        int yP = o1 + (int)(toPixels*(a.r[1]-a.radius));
        int diameterP = (int)(toPixels*a.diameter);
        g.fillOval(xP,yP,diameterP,diameterP);
    }
  }

  //for default boundary widths of zero
  public void initializeCoordinates(double[] simulationRunDimensions) {
    double[] array = {0,0};
    initializeCoordinates(simulationRunDimensions, array);
  }

   // initializeSpecies called by this.paint at design time, and by phase.add at run time
  public void initializeSpecies(Phase phase) {
    setLocation(0,0);
    setSize(phase.getSize());
    double[] simulationRunDimensions;
    if(Beans.isDesignTime()) {
        speciesOrigin[0] = 0;
        speciesOrigin[1] = 0;
        simulationRunDimensions = new double[Space.D];
        simulationRunDimensions[0] = designTimeXDim;
        simulationRunDimensions[1] = designTimeYDim;
    }
    else {
        speciesOrigin[0] = getLocation().x;
        speciesOrigin[1] = getLocation().y;
        simulationRunDimensions = phase.space.getDimensions();
    }
    initializeCoordinates(simulationRunDimensions);
  }

  public void initializeCoordinates(double[] simulationRunDimensions, double[] simulationOrigin) {
    Molecule m = firstMolecule;
    for(Atom a=m.firstAtom; a!=m.lastAtom.getNextAtom(); a=a.getNextAtom()) {
        if(a == m.firstAtom) {a.r[0] = 0.0;}
        else {a.r[0] = a.getPreviousAtom().r[0] + L;}
        a.r[1] = 0.0;
        a.p[0] = a.getMass()*velocity;
        a.p[1] = 0.0;
//        a.velocity = velocity;
    }
    m = m.getNextMolecule();
    for(Atom a=m.firstAtom; a!=m.lastAtom.getNextAtom(); a=a.getNextAtom()) {
        if(a == m.firstAtom) {a.r[0] = 0.0;}
        else {a.r[0] = a.getPreviousAtom().r[0] + L;}
        a.r[1] = simulationRunDimensions[1];
        a.p[0] = -a.getMass()*velocity;
        a.p[1] = 0.0;
    }
  }
}


