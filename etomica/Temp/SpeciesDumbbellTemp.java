package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDumbbell extends Species {

  double velocityX = 10.0, velocityY = 10.0;
  double[] speciesOrigin = {0,0};//origin less depth of boundaries, shifts origin to accomodate boundaries.
  public MoleculeDiatomic[] molecule;
  private boolean showBonds = true;
  double L;
  final double[] e = new double[Space.D];   //default orientation (vertical)
  public static final int nAtoms = 2;
  private transient final int[] pixelDiameter = new int[nAtoms];
  final double[] diameter = new double[nAtoms];

  public SpeciesDumbbell() {
    super();
    diameter[0] = 0.1;
    diameter[1] = 0.1;
    L = 0.10;               //default bondlength
    e[0] = 0.0; e[1] = 1.0; // default orientation (vertical)
    name = "Dumbbell";
    makeMolecules();
  }

  public void makeMolecules() {
    molecule = new MoleculeDiatomic[nElements];
    molecule[0] = new MoleculeDiatomic(speciesIndex,rm,diameter,L,e);
    for(int i = 1; i < nElements; i++) {
        molecule[i] = new MoleculeDiatomic(speciesIndex,rm,diameter,L,e);
        molecule[i-1].setNext(molecule[i]);
        molecule[i].setPrevious(molecule[i-1]);
    }
    firstElement = molecule[0];
    lastElement = molecule[nElements-1];
  }

  // Exposed Properties

  // overrides super classes set method
  public void setNElements(int nElements) {
    this.nElements = nElements;
    makeMolecules();
  }

  public void setDiameter(double d) {   //do away with this when property list can edit arrays
    for (int j=0; j<nAtoms; j++) {this.setDiameter(j,d);}
  }
  public double getDiameter() {return diameter[0];}

  public void setDiameter(int j, double d) {
    this.diameter[j] = d;
    for(int i = 0; i < nElements; i++) {
        molecule[i].setDiameter(j,d);
    }
  }

  public void setSpeciesIndex(int speciesIndex) {
    this.speciesIndex = speciesIndex;
    for(int i = 0; i < nElements; i++) {
        molecule[i].setSpeciesIndex(speciesIndex);
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
  public void setBondLength(double L) {
    this.L = L;
    for(int i=nElements; --i>=0; ) {molecule[i].setBondLength(L);}
  }

  public boolean getShowBonds() {return showBonds;}
  public void setShowBonds(boolean b) {showBonds = b;}

  public void draw(Graphics g, int[] origin, double scale) {
    double toPixels = scale*Phase.TO_PIXELS;
    pixelDiameter[0] = (int)(toPixels*diameter[0]);
    pixelDiameter[1] = (int)(toPixels*diameter[1]);
    g.setColor(color);
    for(int i=nElements; --i>=0; ) {
        for(int j=nAtoms; --j>=0; ) {
          Atom a = molecule[i].atom[j];
          int xP = origin[0] + (int)(toPixels*(a.r[0]-a.getRadius()));
          int yP = origin[1] + (int)(toPixels*(a.r[1]-a.getRadius()));
          g.fillOval(xP,yP,pixelDiameter[j],pixelDiameter[j]);
        }
        if(showBonds) {//draw a line
        }
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
    int moleculeColumns, moleculeRows, i, j, ix, iy;
    double momentumSumX=0, momentumSumY=0, moleculeInitialSpacingX, moleculeInitialSpacingY, momentumNorm;

    //  Number of molecules per row (moleculeColumns) and number of rows (moleculeRows)
    //  in initial configuration
    moleculeColumns = (int)Math.sqrt(((double)simulationRunDimensions[0])/simulationRunDimensions[1]*nElements);
    moleculeRows = (int)(nElements/moleculeColumns);
    //would be better to throw an exception here
    if(moleculeRows*moleculeColumns<nElements)moleculeRows++;
    if(moleculeRows*moleculeColumns < nElements) {
        System.out. println("Program error in SpeciesDisk.initializeCoordinates().moleculeRows");
    }

    //moleculeColumns may be greater than the actual number of columns drawn
    //Need to center columns in the initial position.
    int columnsDrawn = (int)((double)nElements/(double)moleculeRows - 1.0E-10) + 1;
    //moleculeColumnsShift used to center initial coordinates
    double moleculeColumnsShift = simulationRunDimensions[0]/columnsDrawn/2;
    double moleculeRowsShift = simulationRunDimensions[1]/moleculeRows/2;
    //assign distance between molecule centers
    moleculeInitialSpacingX = simulationRunDimensions[0]/columnsDrawn;
    moleculeInitialSpacingY = simulationRunDimensions[1]/moleculeRows;
    Random rand = new Random();
    i = 0;
    double[] rr = new double[Space.D];
    outer: for ( ix=0; ix<moleculeColumns; ix++) {
      inner: for( iy=0; iy<moleculeRows; iy++) {
	    rr[0] = ix*moleculeInitialSpacingX + moleculeColumnsShift + speciesOrigin[0];
	    rr[1] = iy*moleculeInitialSpacingY + moleculeRowsShift + speciesOrigin[1];
        molecule[i].setPosition(rr);
        //assign velocities by random
	    molecule[i].p[1] = Math.cos(2*Math.PI*rand.nextDouble());
	    molecule[i].p[0] = Math.sqrt(1.0 - molecule[i].p[1]*molecule[i].p[1]);

	    momentumNorm = Math.sqrt(molecule[i].p[0]*molecule[i].p[0]+molecule[i].p[1]*molecule[i].p[1]);
	    molecule[i].p[0] *= velocityX/momentumNorm/rm/sigma;
	    molecule[i].p[1] *= velocityY/momentumNorm/rm/sigma;
	    momentumSumX += molecule[i].p[0]; momentumSumY += molecule[i].p[1];
	    i++;
	if(i == nElements) {break outer;}
      }
    }

//    Zero center-of-mass momentum
    momentumSumX /= nElements; momentumSumY /= nElements;
    double[] e = new double[Space.D];
    for(i=0; i<nElements; i++) {
        molecule[i].p[0] -= momentumSumX;
        molecule[i].p[1] -= momentumSumY;

//        molecule[i].r[0] += 0.01*(2.0*rand.nextDouble()-1.0);
//        molecule[i].r[1] += 0.01*(2.0*rand.nextDouble()-1.0);
        }
  }
}


