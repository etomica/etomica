package simulate;
import java.util.*;
import java.awt.*;

public class Configuration1 extends Configuration{

    private double temperature = 300;
    
    public void initializeCoordinates() {
        int sumOfMolecules = 0;
        Species s;
        for(Enumeration e = species.elements(); e.hasMoreElements(); ){
            s = (Species)e.nextElement();
            sumOfMolecules+=s.getNMolecules();
        }
        s = (Species)species.firstElement();
        double[] simulationOrigin = {0,0};
        int moleculeColumns, moleculeRows, i, j, ix, iy;
        double momentumSumX=0, momentumSumY=0, moleculeInitialSpacingX, moleculeInitialSpacingY, momentumNorm;
        double momentum = Math.sqrt(s.mass*temperature/Constants.KE2T*(double)Space.D);  //need to divide by sqrt(m) to get velocity
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
        ix = iy = 0;
        for(Enumeration ee = species.elements(); ee.hasMoreElements(); ){
            s = (Species)ee.nextElement();
            Molecule m = s.firstMolecule;
            outer: for ( ; ix<moleculeColumns;ix++) {
            if (iy >= moleculeRows) {iy=0;}
            inner: for( ; iy<moleculeRows;iy++) {
                Atom a = m.atom[0];
                a.r[0] = ix*moleculeInitialSpacingX + moleculeColumnsShift + 0.5*s.nAtomsPerMolecule*e[0]*s.L;
	            a.r[1] = iy*moleculeInitialSpacingY + moleculeRowsShift + 0.5*s.nAtomsPerMolecule*e[1]*s.L;
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
                if(m == s.lastMolecule) {
                    if(iy+1==moleculeRows){
                       ix++;}
                    iy++;
                    break outer;}
                m=m.getNextMolecule();
                }
            }
        }

    //    Zero center-of-mass momentum
        momentumSumX /= sumOfMolecules; momentumSumY /= sumOfMolecules;
        for(i=0; i<sumOfMolecules; i++) {
//            m.atom[0].p[0] -= momentumSumX;
//            m.atom[0].p[1] -= momentumSumY;
            }
      }

    public double getTemperature(){
        return temperature;
    }
    
    public void setTemperature(double t){
        temperature = t;
        initializeCoordinates();
    }
}
