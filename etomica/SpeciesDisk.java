package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDisk extends Species {

  double[] speciesOrigin = {0,0}; //origin less depth of boundaries, shifts origin to accomodate boundaries.
  private transient final int[] shiftOrigin = new int[Space.D];     //work vector for drawing overflow images
  double mass;
  double diameter;
  double radius;
  Color color;
  
  public SpeciesDisk() {
    super();
  }
  
  public void setDefaults() {
    setDiameter(0.1);    
    setMass(1.0);
    setColor(Color.black);
    name = "Disk";
    nAtomsPerMolecule = 1;
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

    public final double getMass() {return mass;}
    public final void setMass(double mass) {
        this.mass = mass;
        if(firstAtom == null) {return;}  //return if atoms have not yet been ordered
        for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {a.setMass(mass);}
        for(Molecule m=firstMolecule; m!=lastMolecule.getNextMolecule(); m=m.getNextMolecule()) {m.updateMass();}        
    }
    
    public final double getDiameter() {return diameter;}
    public void setDiameter(double d) {
        diameter = d;
        radius = 0.5*d;
        if(firstAtom == null) {return;}
        for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {a.setDiameter(d);}
    }
        
    public final Color getColor() {return color;}
    public final void setColor(Color c) {
        color = c;
        if(firstAtom == null) {return;}
        for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {a.setColor(c);}
    }
        
  public void draw(Graphics g, int[] origin, double scale) {

    double toPixels = scale*Phase.TO_PIXELS;
 /*   
    int diameterP = (int)(toPixels*diameter);
    g.setColor(color);
    */
    Atom nextSpeciesAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
        a.draw(g,origin,scale);
        /*
        int xP = origin[0] + (int)(toPixels*(a.r[0]-radius));
        int yP = origin[1] + (int)(toPixels*(a.r[1]-radius));
        g.fillOval(xP,yP,diameterP,diameterP);
        */
    }
    if(parentPhase.drawOverflowImages) {
        for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
            double[][] shifts = parentPhase.space.getOverflowShifts(a.r,radius);
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
/*    if(parentPhase.useNeighborList) {
        for(Molecule m=firstMolecule; m!=lastMolecule.getNextMolecule(); m=m.getNextMolecule()) {
            m.needNeighborUpdate = true;
        }
    }
    */
  }

  public void initializeCoordinates(double[] simulationRunDimensions, double[] simulationOrigin) {
    int moleculeColumns, moleculeRows, i, j, ix, iy;
    double momentumSumX=0, momentumSumY=0, moleculeInitialSpacingX, moleculeInitialSpacingY, momentumNorm;
    double momentum = Math.sqrt(mass*parentPhase.getInitialTemperature()/Constants.KE2T*(double)Space.D);  //need to divide by sqrt(m) to get velocity
    //  Number of molecules per row (moleculeColumns) and number of rows (moleculeRows)
    //  in initial configuration
    moleculeColumns = (int)Math.sqrt(((double)simulationRunDimensions[0])/simulationRunDimensions[1]*nMolecules);
    moleculeRows = (int)(nMolecules/moleculeColumns);
    //would be better to throw an exception here
    if(moleculeRows*moleculeColumns<nMolecules)moleculeRows++;
    if(moleculeRows*moleculeColumns < nMolecules) {
        System.out.println("Program error in SpeciesDisk.initializeCoordinates().moleculeRows");
    }
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
	    a.r[0] = ix*moleculeInitialSpacingX + moleculeColumnsShift + speciesOrigin[0];
	    a.r[1] = iy*moleculeInitialSpacingY + moleculeRowsShift + speciesOrigin[1];
        //assign velocities by random
	    a.p[1] = Math.cos(2*Math.PI*rand.nextDouble());
	    a.p[0] = Math.sqrt(1.0 - molecule[i].p[1]*molecule[i].p[1]);

	    //debugging code
/*	    molecule[i].p[0] = 0.0;
	    if(molecule[i].r[1] > 0.5) {molecule[i].p[1] = -1.0;}
	    else {molecule[i].p[1] = +1.0;}
*/	    //debugging code
	    
	    momentumNorm = Math.sqrt(a.p[0]*a.p[0]+a.p[1]*a.p[1]);
	    a.p[0] *= momentum/momentumNorm; 
	    a.p[1] *= momentum/momentumNorm;
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


