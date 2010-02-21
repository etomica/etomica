package etomica.models.oneDHardRods;

import etomica.api.IVectorMutable;
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

public class DoubleIntegralEtas extends Simulation{
    
    NormalModes nm;
    protected CoordinateDefinition cDef;
    private double[][][] eigenvectors;
    private IVectorMutable[] waveVectors;
    private double[] wvc;
    double[][] omega;

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
        
        cDef = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        cDef.initializeCoordinates(nCells);
        
        bCells = cDef.getBasisCells();
        coordinateDim = cDef.getCoordinateDim();
        normalization = 1/Math.sqrt(bCells.length);
        sqrtT = Math.sqrt(1.0);
        atomLocs = new double[bCells.length];
        u = new double[bCells.length];
        x0Pos = new double[bCells.length];
        x0Pos[0] = -2.142857142857143;
        x0Pos[1] = -0.7142857142857142;
        x0Pos[2] = 0.7142857142857142;


        nm = new NormalModes1DHR(box.getBoundary(), nAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(1.0);
        nm.getOmegaSquared();
        nm.getWaveVectorFactory().makeWaveVectors(box);
        waveVectors = nm.getWaveVectorFactory().getWaveVectors();
        
        wvc = nm.getWaveVectorFactory().getCoefficients();
        omega = nm.getOmegaSquared();
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
        if(Math.abs(function(xValue, yValue)) >= tol){
            System.out.println("Increase your ranges xStart, yStart.");
        }
            
        yValue = yEnd;
        if(Math.abs(function(xValue, yValue)) >= tol){
            System.out.println("Increase your ranges xStart, yEnd.");
        }
        xValue = xEnd;
        yValue = yStart;
        if(Math.abs(function(xValue, yValue)) >= tol){
            System.out.println("Increase your ranges xEnd, yStart.");
        }
        
        yValue = yEnd;
        if(Math.abs(function(xValue, yValue)) >= tol){
            System.out.println("Increase your ranges xEnd, yEnd.");
        }
        
        //Now do the "edge" stuff, where either the x or y value, but not both,
        //  are the start or end value
        xValue = xStart;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            if(Math.abs(function(xValue, yValue)) >= tol){
                System.out.println("Increase your ranges xStart, yEdge.");
            }
        }
        xValue = xEnd;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            total += 2 * function(xValue, yValue);
            if(Math.abs(function(xValue, yValue)) >= tol){
                System.out.println("Increase your ranges xEnd, yEdge.");
            }
        }
        yValue = yStart;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;if(Math.abs(function(xValue, yValue)) >= tol){
                System.out.println("Increase your ranges xEdge, yStart.");
            }
        }
        yValue = yEnd;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(function(xValue, yValue)) >= tol){
                System.out.println("Increase your ranges xEdge, yEnd.");
            }
        }
        
        
        //The actual calculation of the integral.
        //  We have already eliminated the edges' and endpoints' contributions,
        //  because they have been confirmed to be zero by the above code.
        for(int i = 1; i < xN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            for (int j = 1; j < yN; j++) {
                yValue = yStart + j *(yEnd - yStart) / yN;
                if(!overlap(xValue, yValue)){
                    total += 4 * function(xValue, yValue);
                }
            }
        }
        
        //Now we do the prefix thing
        double prefix = (xEnd - xStart)*(yEnd - yStart) / (4 * xN * yN);
        total *= prefix;
        
        System.out.println("Integral = " + total);
        return total;
    }
    
    private boolean overlap(double etaReal, double etaImag){
        boolean isOverlap = false;
        
        //This pile o' stuff here calculates the x position of each atom.
        for(int iCell = 0; iCell < bCells.length; iCell++ ){
            BasisCell cell = bCells[iCell];
            for(int i = 0; i < coordinateDim; i++){
                u[i] = 0.0;
            }
            
            for(int iVector = 0; iVector < waveVectors.length; iVector++){
                double kR = waveVectors[iVector].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                
                for (int i=0; i<coordinateDim; i++) {
                    for (int j=0; j<coordinateDim; j++) {
                        u[j] += wvc[iVector] * eigenvectors[iVector][i][j] * 2.0 *
                            (etaReal * coskR + etaImag * sinkR);
                        
                    }
                }
            }
            
            for(int i=0; i < coordinateDim; i++){
                u[i] *= normalization;
            }
            
            atomLocs[iCell] = u[iCell] + x0Pos[iCell];
        }
        
        for (int i=0; i < atomLocs.length-1; i++){
            if ( (atomLocs[i] + 0.5) >= (atomLocs[i+1] -0.5) ){
                isOverlap = true;
                break;
            }
        }
        
        double repeat = box.getBoundary().getEdgeVector(0).getX(0);
        if( (atomLocs[0]+repeat-0.5) <= (atomLocs[atomLocs.length-1]+0.5) ){
            isOverlap = true;
        }
        
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

    private double function(double x, double y){
        double value = x * y ;
        return value;
    }
    
    public static void main(String[] args) {
      double xStart = -1.0;
      double yStart = -1.0;
      double xEnd = 1.0;
      double yEnd = 1.0;
      int xN = 100000;
      int yN = 100000;
      
      int nAtoms = 3;
      double density = 0.7;
      
      DoubleIntegralEtas di = new DoubleIntegralEtas(nAtoms, density);
      di.setIntegrationParameters(xStart, xEnd, yStart, yEnd, xN, yN);
      di.calculate();
    }
}
