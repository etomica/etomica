package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDumbbell extends Species {

  double velocityX = 50.0, velocityY = 50.0;
  double[] speciesOrigin = {0.,0.}; //origin less depth of boundaries, shifts origin to accomodate boundaries.
  double[] mass;      //make arrays when doublearrayeditor in place
  double[] diameter;
  double[] radius;
  Color color;
  
  private boolean showBonds = true;
  double L;    //atom-atom separation
  final double[] e = new double[Space.D];   //default orientation (vertical)
  public static final int nAtomsPerMolecule = 2;
 
  public SpeciesDumbbell() {
    super();
    setColor(Color.black);
    L = 0.10;               //default bondlength
    e[0] = 0.0; e[1] = 1.0; // default orientation (vertical)
    name = "Dumbbell";
  }
  
  public void setDefaults() {
    mass = new double[nAtomsPerMolecule];
    diameter = new double[nAtomsPerMolecule];
    radius = new double[nAtomsPerMolecule];
    for(int i=0; i<nAtomsPerMolecule; i++) {
        mass[i] = 1.0;
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
 
  public double getVelocityX() {return velocityX;}
  public void setVelocityX(double x){
    velocityX = x;
  }

  public double getVelocityY() {return velocityY;}
  public void setVelocityY(double y){
    velocityY = y;
  }

  public double getBondLength() {return L;}
  public void setBondLength(double L) {this.L = L;}

  public boolean getShowBonds() {return showBonds;}
  public void setShowBonds(boolean b) {showBonds = b;}

  public void draw(Graphics g, int[] origin, double scale) {
    double toPixels = scale*Phase.TO_PIXELS;
    g.setColor(color);
    int o0 = origin[0];  //copy to scalars to avoid array-bounds check in loop
    int o1 = origin[1];
    for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {
        int xP = o0 + (int)(toPixels*(a.r[0]-a.radius));
        int yP = o1 + (int)(toPixels*(a.r[1]-a.radius));
        int diameterP = (int)(toPixels*a.diameter);
        g.fillOval(xP,yP,diameterP,diameterP);
    }
    if(showBonds) {//loop over molecules and draw a line
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
        speciesOrigin[0] = 0;
        speciesOrigin[1] = 0;
        simulationRunDimensions = new double[Space.D];
        simulationRunDimensions[0] = designTimeXDim;
        simulationRunDimensions[1] = designTimeYDim;
    }
    else {
        speciesOrigin[0] = getLocation().x;
        speciesOrigin[1] = getLocation().y;
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
//        for(Molecule m=firstMolecule; m!=lastMolecule.getNextMolecule(); m=m.getNextMolecule()) {
//            m.needNeighborUpdate = true;
//        }
  }

  public void initializeCoordinates(double[] simulationRunDimensions, double[] simulationOrigin) {
    int moleculeColumns, moleculeRows, i, j, ix, iy;
    double momentumSumX=0, momentumSumY=0, moleculeInitialSpacingX, moleculeInitialSpacingY, momentumNorm;

    //  Number of molecules per row (moleculeColumns) and number of rows (moleculeRows)
    //  in initial configuration
    moleculeColumns = (int)Math.sqrt(((double)simulationRunDimensions[0])/simulationRunDimensions[1]*nMolecules);
    moleculeRows = (int)(nMolecules/moleculeColumns);
    //would be better to throw an exception here
    if(moleculeRows*moleculeColumns<nMolecules)moleculeRows++;
    if(moleculeRows*moleculeColumns < nMolecules) {
        System.out.println("Program error in SpeciesDisk.initializeCoordinates().moleculeRows");
    }
//   double sig = numberDensity*Math.min(simulationRunDimensions[0]/(double)moleculeColumns,simulationRunDimensions[1]/(double)moleculeRows);
//   setSigma(sig);

    //moleculeColumns may be greater than the actual number of columns drawn
    //Need to center columns in the initial position.
    int columnsDrawn = (int)((double)nMolecules/(double)moleculeRows - 1.0E-10) + 1;
    //moleculeColumnsShift used to center initial coordinates
    double moleculeColumnsShift = simulationRunDimensions[0]/columnsDrawn/2;    
    double moleculeRowsShift = simulationRunDimensions[1]/moleculeRows/2;
    //assign distance between molecule centers 
    moleculeInitialSpacingX = simulationRunDimensions[0]/columnsDrawn;
    moleculeInitialSpacingY = simulationRunDimensions[1]/moleculeRows;
    Random rand = new Random();    
    i = 0;
    outer: for ( ix=0; ix<moleculeColumns; ix++) {
      inner: for( iy=0; iy<moleculeRows; iy++) {
        Atom a = molecule[i].atom[0];
	    a.r[0] = ix*moleculeInitialSpacingX + moleculeColumnsShift + speciesOrigin[0] + 0.5*e[0]*L;
	    a.r[1] = iy*moleculeInitialSpacingY + moleculeRowsShift + speciesOrigin[1] + 0.5*e[1]*L;
        a = molecule[i].atom[1];
	    a.r[0] = ix*moleculeInitialSpacingX + moleculeColumnsShift + speciesOrigin[0] - 0.5*e[0]*L;
	    a.r[1] = iy*moleculeInitialSpacingY + moleculeRowsShift + speciesOrigin[1] - 0.5*e[1]*L;
        //assign velocities by random
	    a.p[1] = Math.cos(2*Math.PI*rand.nextDouble());
	    a.p[0] = Math.sqrt(1.0 - molecule[i].p[1]*molecule[i].p[1]);

	    momentumNorm = Math.sqrt(a.p[0]*a.p[0]+a.p[1]*a.p[1]);
	    a.p[0] *= mass[1]*velocityX/momentumNorm/a.diameter; 
	    a.p[1] *= mass[1]*velocityY/momentumNorm/a.diameter;
	    momentumSumX += a.p[0]; momentumSumY += a.p[1];
	    i++;
	if(i == nMolecules) {break outer;}
      }
    }

//    Zero center-of-mass momentum
    momentumSumX /= nMolecules; momentumSumY /= nMolecules;
    for(i=0; i<nMolecules; i++) {
        molecule[i].atom[0].p[0] -= momentumSumX; 
        molecule[i].atom[0].p[1] -= momentumSumY;
//        molecule[i].r[0] += 0.01*(2.0*rand.nextDouble()-1.0);
//        molecule[i].r[1] += 0.01*(2.0*rand.nextDouble()-1.0);
        }
  }
}


