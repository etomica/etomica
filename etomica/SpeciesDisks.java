package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDisks extends Species {
 ///     public double L;    //atom-atom separation, used only for initial configuration;


  double[] speciesOrigin = {0,0}; //origin less depth of boundaries, shifts origin to accomodate boundaries.
///  double mass;
  double diameter;
  double radius;
  Color color;
  
  public SpeciesDisks() {
    super();
    initializeMolecules(0.1,1.0,Color.black);
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

    public void setL(double L) {this.L = L;}
    public double getL() {return L;}

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
        

  //for default boundary widths of zero
 /// public void initializeCoordinates(double[] simulationRunDimensions) {
 ///   double[] array = {0,0};
 ///   initializeCoordinates(simulationRunDimensions, array);
 /// }

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
 ///   initializeCoordinates(simulationRunDimensions);
/*    if(parentPhase.useNeighborList) {
        for(Molecule m=firstMolecule; m!=lastMolecule.getNextMolecule(); m=m.getNextMolecule()) {
            m.needNeighborUpdate = true;
        }
    }
    */
  }
/*
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
    double e[] = new double[2];
    e[0] = 0.71;
    e[1] = 0.71;
    outer: for ( ix=0; ix<moleculeColumns; ix++) {
      inner: for( iy=0; iy<moleculeRows; iy++) {
        Atom a = molecule[i].atom[0];
	    a.r[0] = ix*moleculeInitialSpacingX + moleculeColumnsShift + speciesOrigin[0] + 0.5*nAtomsPerMolecule*e[0]*L;
	    a.r[1] = iy*moleculeInitialSpacingY + moleculeRowsShift + speciesOrigin[1] + 0.5*nAtomsPerMolecule*e[1]*L;
	    for(int jj=1; jj<nAtomsPerMolecule; jj++) {
            a = molecule[i].atom[jj];
	        a.r[0] = a.getPreviousAtom().r[0] - e[0]*L;
	        a.r[1] = a.getPreviousAtom().r[1] - e[1]*L;
	    }
        //assign velocities by random
	    a.p[1] = Math.cos(2*Math.PI*rand.nextDouble());
	    a.p[0] = Math.sqrt(1.0 - molecule[i].p[1]*molecule[i].p[1]);

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
  }*/
}


