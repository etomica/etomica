package simulate;
import java.util.*;
import java.awt.*;

public class Configuration2 extends Configuration{

    public void initializeCoordinates(){
        int sumOfMolecules = 0;
        int i = 0;
        Species s;
        for(Enumeration e = species.elements(); e.hasMoreElements(); ){
            s = (Species)e.nextElement();
            sumOfMolecules+=s.getNMolecules();
        }
        for(Enumeration e = species.elements(); e.hasMoreElements(); ){
            s = (Species)e.nextElement();
            initializeCoordinates(i,sumOfMolecules,s);
            i++;
        }
        s = (Species)species.firstElement();
        s.parentPhase.repaint();   
    }
    public void initializeCoordinates(int index, int sumOfMolecules, Species s) {
        double[] simulationOrigin = {0,0};
        int moleculeColumns, moleculeRows, i, j, ix, iy;
        double momentumSumX=0, momentumSumY=0, moleculeInitialSpacingX, moleculeInitialSpacingY, momentumNorm;
        double momentum = Math.sqrt(s.mass*s.parentPhase.getInitialTemperature()/Constants.KE2T*(double)Space.D);  //need to divide by sqrt(m) to get velocity
    //  Number of molecules per row (moleculeColumns) and number of rows (moleculeRows)
    //  in initial configuration
        moleculeColumns = (int)Math.sqrt(((double)s.simulationRunDimensions[0])/s.simulationRunDimensions[1]*sumOfMolecules);
        moleculeRows = (int)(sumOfMolecules/moleculeColumns);
    //would be better to throw an exception here
        if(moleculeRows*moleculeColumns<sumOfMolecules)moleculeRows++;
        if(moleculeRows*moleculeColumns < sumOfMolecules) {
            System.out.println("Program error in SpeciesDisk.initializeCoordinates().moleculeRows");
        }
    //moleculeColumns may be greater than the actual number of columns drawn
    //Need to center columns in the initial position.
        int columnsDrawn = (int)((double)sumOfMolecules/(double)moleculeRows - 1.0E-10) + 1;
    //moleculeColumnsShift used to center initial coordinates
        double moleculeColumnsShift = s.simulationRunDimensions[0]/columnsDrawn/2;
        double moleculeRowsShift = s.simulationRunDimensions[1]/moleculeRows/2;
    //assign distance between molecule centers
        moleculeInitialSpacingX = s.simulationRunDimensions[0]/columnsDrawn;
        moleculeInitialSpacingY = s.simulationRunDimensions[1]/moleculeRows;
        Random rand = new Random();
        i = 0;
        double e[] = new double[2];
        e[0] = 0.71;
        e[1] = 0.71;
        Molecule m = s.molecule[i];
        outer: for ( ix=0; ix<moleculeColumns; ix++) {
          inner: for( iy=0; iy<moleculeRows; iy++) {
            if((ix+iy)>0)
                m=m.getNextMolecule();
            Atom a = m.atom[0];
            a.r[0] = ix*moleculeInitialSpacingX + (2-(double)s.getNMolecules()/sumOfMolecules)*(ix+index+1)*moleculeColumnsShift + 0.5*s.nAtomsPerMolecule*e[0]*s.L;
	        a.r[1] = iy*moleculeInitialSpacingY + (2-(double)s.getNMolecules()/sumOfMolecules)*(index+1)*moleculeRowsShift + 0.5*s.nAtomsPerMolecule*e[1]*s.L;
    	    for(int jj=1; jj<s.nAtomsPerMolecule; jj++) {
                a = m.atom[jj];
	            a.r[0] = a.getPreviousAtom().r[0] - e[0]*s.L;
	            a.r[1] = a.getPreviousAtom().r[1] - e[1]*s.L;
	        }
        //assign velocities by random
	        a.p[1] = Math.cos(2*Math.PI*rand.nextDouble());
    	    a.p[0] = Math.sqrt(1.0 - m.p[1]*m.p[1]);

	        momentumNorm = Math.sqrt(a.p[0]*a.p[0]+a.p[1]*a.p[1]);
	        a.p[0] *= momentum/momentumNorm;
    	    a.p[1] *= momentum/momentumNorm;
	        momentumSumX += a.p[0]; momentumSumY += a.p[1];
	        i++;
    	if(i == s.getNMolecules()) {break outer;}
          }
        }

    //    Zero center-of-mass momentum
        momentumSumX /= s.getNMolecules(); momentumSumY /= s.getNMolecules();
        for(i=0; i<s.getNMolecules(); i++) {
            s.molecule[i].atom[0].p[0] -= momentumSumX;
            s.molecule[i].atom[0].p[1] -= momentumSumY;
        }
    }
    
}
