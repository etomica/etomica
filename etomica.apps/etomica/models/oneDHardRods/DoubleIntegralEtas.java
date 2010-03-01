package etomica.models.oneDHardRods;

import etomica.api.IAtomList;
import etomica.api.IVectorMutable;
import etomica.atom.Atom;
import etomica.box.Box;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.NormalModes;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * Please note this class is specifically designed to integrate for 2 or 3
 * rod systems, with one rod per cell.
 * 
 * Hardcoded with the initial atom locations for 2 or 3 atoms.
 * 
 * Hardcoded to only use the wave vector #1
 * 
 * @author cribbin
 *
 */
public class DoubleIntegralEtas extends Simulation{
    
    NormalModes nm;
    protected CoordinateDefinition cDef;
    private double[][][] eigenvectors;
    private IVectorMutable[] waveVectors;
    private double[] wvc;
    double[][] omega2;

    Box box;
    double xStart, xEnd, yStart, yEnd;
    int xN, yN, nAtoms;
    double density;
    int[] nCells;
    
    BasisCell[] bCells;
    int coordinateDim;
    double normalization;
    double sqrtT;
    double[] atomLocs, u, x0Pos;

    
    public DoubleIntegralEtas(int nAtoms, double density){
        super(Space.getInstance(1));
        
        this.nAtoms = nAtoms;
        this.density = density;
        
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        Basis basis = new BasisMonatomic(space);
        
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, nAtoms);
        
        Primitive primitive = new PrimitiveCubic(space, 1.0/density);
        Boundary bdry = new BoundaryRectangularPeriodic(space, nAtoms/density);
        nCells = new int[]{nAtoms};
        box.setBoundary(bdry);
        
        cDef = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        cDef.initializeCoordinates(nCells);
        
        bCells = cDef.getBasisCells();
        coordinateDim = cDef.getCoordinateDim();
        normalization = 1/Math.sqrt(bCells.length);
        sqrtT = Math.sqrt(1.0);
        atomLocs = new double[bCells.length];
        u = new double[bCells.length];
        x0Pos = new double[bCells.length];

        //THESE ARE HARD CODE FOR 2 RODS, 0.7 density
        if(nAtoms == 2){
            x0Pos[0] = -1.4285714285714286;
            x0Pos[1] = 0.0;
        }
        if(nAtoms == 3){
            x0Pos[0] = -2.142857142857143;
            x0Pos[1] = -0.7142857142857142;
            x0Pos[2] = 0.7142857142857142;
        }
            
        nm = new NormalModes1DHR(box.getBoundary(), nAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(1.0);
        nm.getOmegaSquared();
        nm.getWaveVectorFactory().makeWaveVectors(box);
        waveVectors = nm.getWaveVectorFactory().getWaveVectors();
        
        wvc = nm.getWaveVectorFactory().getCoefficients();
        omega2 = nm.getOmegaSquared();
        eigenvectors = nm.getEigenvectors();
    }
    
    
    public double calculate(){
        double total = 0.0;
        double xValue = 0.0;
        double yValue = 0.0;
        
        //Here we are checking that we have enclosed all values of potential
        //  interest.
        //First, the four values for which both the x and y values are
        //  the start or end value
        double tol = 0.0000000001;
        xValue = xStart;
        yValue = yStart;
        if(Math.abs(integrand(xValue, yValue)) >= tol ){
            System.out.println("Increase your ranges xStart, yStart.");
        }
            
        yValue = yEnd;
        if(Math.abs(integrand(xValue, yValue)) >= tol ){
            System.out.println("Increase your ranges xStart, yEnd.");
        }
        xValue = xEnd;
        yValue = yStart;
        if(Math.abs(integrand(xValue, yValue)) >= tol ){
            System.out.println("Increase your ranges xEnd, yStart.");
        }
        
        yValue = yEnd;
        if(Math.abs(integrand(xValue, yValue)) >= tol ){
            System.out.println("Increase your ranges xEnd, yEnd.");
        }
        
        //Now do the "edge" stuff, where either the x or y value, but not both,
        //  are the start or end value
        xValue = xStart;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            if(Math.abs(integrand(xValue, yValue)) >= tol ){
                System.out.println("Increase your ranges xStart, yEdge.");
            }
        }
        xValue = xEnd;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            total += 2 * integrand(xValue, yValue);
            if(Math.abs(integrand(xValue, yValue)) >= tol ){
                System.out.println("Increase your ranges xEnd, yEdge.");
            }
        }
        yValue = yStart;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(integrand(xValue, yValue)) >= tol ){
                        System.out.println("Increase your ranges xEdge, yStart. " + integrand(xValue, yValue));
            }
        }
        yValue = yEnd;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(integrand(xValue, yValue)) >= tol ){
                System.out.println("Increase your ranges xEdge, yEnd " + integrand(xValue, yValue));
            }
        }
        
        
        //The actual calculation of the integral.
        //  We have already eliminated the edges' and endpoints' contributions,
        //  because they have been confirmed to be zero by the above code.
        for(int i = 1; i < xN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            for (int j = 1; j < yN; j++) {
                yValue = yStart + j *(yEnd - yStart) / yN;
                    total += 4 * integrand(xValue, yValue);
            }
        }
        //nan could also code the above with xValue = xStart; xValue  += (xEnd - xStart) / xN;
        
        //Now we do the prefix thing
        double prefix = (xEnd - xStart)*(yEnd - yStart) / (4 * xN * yN);
        total *= prefix;
        System.out.println("Integral = " + total);
        return total;
    }
    
    private boolean overlap(double etaReal, double etaImag){
        boolean isOverlap = false;
        
        //This pile o' stuff here calculates the x position of each atom.
        for(int i = 0; i < nAtoms; i++){
            u[i] = 0.0;
        }
        
        for(int iCell = 0; iCell < bCells.length; iCell++ ){
            BasisCell cell = bCells[iCell];
            
            double kR = waveVectors[1].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);

            if( !(Double.isInfinite(omega2[1][0])) ){
                if(nAtoms == 2) {
                    u[iCell] += wvc[1] * eigenvectors[1][0][0] * 2.0 * (etaReal * coskR);
                }else{
                    u[iCell] += wvc[1] * eigenvectors[1][0][0] * 2.0 * (etaReal * coskR - etaImag * sinkR);
                }
            }
            
            u[iCell] *= normalization;
            atomLocs[iCell] = u[iCell] + x0Pos[iCell];
            
        }//end of iCell loop
        
        for (int i=0; i < atomLocs.length-1; i++){
            if ( (atomLocs[i] + 0.5) >= (atomLocs[i+1] - 0.5) ){
                isOverlap = true;
                break;
            }
        }
        
        double repeat = box.getBoundary().getBoxSize().getX(0);
        if( (atomLocs[0]+repeat-0.5) <= (atomLocs[atomLocs.length-1]+0.5) ){
            isOverlap = true;
        }
        
//        if (!isOverlap) {System.out.println(etaReal + " " + etaImag);}
        System.out.println("etaReal " + etaReal + " " + isOverlap + " 0 " + atomLocs[0]+" 1 " + atomLocs[1]);
       
        return isOverlap;
    }
    
    
    public void setIntegrationParameters(double xStart, double xEnd, double yStart, double yEnd,
            int xN, int yN){
        this.xStart = xStart;
        this.xEnd = xEnd;
        this.yStart = yStart;
        this.yEnd = yEnd;
        this.xN = xN;
        this.yN = yN;
    }

    private double integrand(double x, double y){
        double value;
        if (overlap(x, y)) {
            value = 0;
        }else{
            value = 1;
        }
        return value;
    }
    
    public static void main(String[] args) {
      double xStart = -1.0;
      double yStart = -1.0;
      double xEnd = 1.0;
      double yEnd = 1.0;
      
      int xN = 100;
      int yN = 2;
      
      int nAtoms = 2;
      double density = 0.7;
      
      DoubleIntegralEtas di = new DoubleIntegralEtas(nAtoms, density);
      di.setIntegrationParameters(xStart, xEnd, yStart, yEnd, xN, yN);
      di.calculate();
    }
}
