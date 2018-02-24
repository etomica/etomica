/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.space.Vector;
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
    
    double sqrtT;
    double density;
    Box boxT, boxR;
    double alpha;
    double xStart, xEnd, yStart, yEnd;
    int xN, yN, nAtomsT, nAtomsR;
    NormalModes nmT, nmR;
    protected CoordinateDefinition cDefT, cDefR;
    private Vector[] wvT, wvR;
    private double[] wvcT, sqrtWvcT, wvcR, sqrtWvcR;
    double[][] omega2T, oneOverOmega2T, omega2R, oneOverOmega2R;
    int[] nCellsT, nCellsR;
    BasisCell[] rCells, tCells;
    double normalizationT, normalizationR;
    double[] atomLocsT, uT, x0PosT, atomLocsR, uR, x0PosR;

    
    public DoubleIntegralEtas(int nAtoms, double density, double a) {
        super(Space.getInstance(1));
        this.density = density;
        sqrtT = Math.sqrt(1.0);
        alpha = a;

        //TARGET
        this.nAtomsT = nAtoms;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        Basis basis = new BasisMonatomic(space);

        boxT = this.makeBox();
        boxT.setNMolecules(species, nAtomsT);

        Primitive primitive = new PrimitiveCubic(space, 1.0 / density);
        Boundary bdryT = new BoundaryRectangularPeriodic(space, nAtomsT / density);
        nCellsT = new int[]{nAtomsT};
        boxT.setBoundary(bdryT);

        cDefT = new CoordinateDefinitionLeaf(boxT, primitive, basis, space);
        cDefT.initializeCoordinates(nCellsT);

        tCells = cDefT.getBasisCells();
        normalizationT = 1 / Math.sqrt(tCells.length);
        atomLocsT = new double[tCells.length];
        uT = new double[tCells.length];
        x0PosT = new double[tCells.length];

        x0PosT[0] = -2.142857142857143;
        x0PosT[1] = -0.7142857142857142;
        x0PosT[2] = 0.7142857142857142;

        nmT = new NormalModes1DHR(boxT.getBoundary(), nAtoms);
        nmT.setHarmonicFudge(1.0);
        nmT.setTemperature(1.0);
        nmT.getOmegaSquared();
        nmT.getWaveVectorFactory().makeWaveVectors(boxT);
        wvT = nmT.getWaveVectorFactory().getWaveVectors();

        wvcT = nmT.getWaveVectorFactory().getCoefficients();
        sqrtWvcT = new double[wvcT.length];
        for (int i = 0; i < wvcT.length; i++) {
            sqrtWvcT[i] = Math.sqrt(2 * wvcT[i]);
        }
        omega2T = nmT.getOmegaSquared();
        oneOverOmega2T = new double[omega2T.length][omega2T[0].length];
        for (int i = 0; i < omega2T.length; i++) {
            for (int j = 0; j < omega2T[i].length; j++) {
                oneOverOmega2T[i][j] = Math.sqrt(1.0 / (omega2T[i][j]));
            }
        }

        //REFERENCE
        nAtomsR = nAtomsT - 1;
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        basis = new BasisMonatomic(space);

        Boundary bdryR = new BoundaryRectangularPeriodic(space, nAtomsR / density);
        boxR = this.makeBox(bdryR);
        boxR.setNMolecules(species, nAtomsR);

        primitive = new PrimitiveCubic(space, 1.0 / density);
        nCellsR = new int[]{nAtomsR};

        cDefR = new CoordinateDefinitionLeaf(boxR, primitive, basis, space);
        cDefR.initializeCoordinates(nCellsR);

        rCells = cDefR.getBasisCells();
        normalizationR = 1 / Math.sqrt(rCells.length);
        atomLocsR = new double[rCells.length];
        uR = new double[rCells.length];
        x0PosR = new double[rCells.length];

        //THESE ARE HARD CODE FOR 2 RODS, 0.7 density
        x0PosR[0] = -1.4285714285714286;
        x0PosR[1] = 0.0;

        nmR = new NormalModes1DHR(boxR.getBoundary(), nAtomsR);
        nmR.setHarmonicFudge(1.0);
        nmR.setTemperature(1.0);
        nmR.getWaveVectorFactory().makeWaveVectors(boxR);
        wvR = nmR.getWaveVectorFactory().getWaveVectors();

        wvcR = nmR.getWaveVectorFactory().getCoefficients();
        sqrtWvcR = new double[wvcR.length];
        for (int i = 0; i < wvcR.length; i++) {
            sqrtWvcR[i] = Math.sqrt(2 * wvcR[i]);
        }
        omega2R = nmT.getOmegaSquared();
        oneOverOmega2R = new double[omega2R.length][omega2R[0].length];
        for (int i = 0; i < omega2R.length; i++) {
            for (int j = 0; j < omega2R[i].length; j++) {
                oneOverOmega2R[i][j] = Math.sqrt(1.0 / (omega2R[i][j]));
            }
        }


    }
    
    public double calculate(){
        double numerator = getNumerator();
        double denomRef = getDenomRef();
        double denomTarg = getDenomTarg();
        
        double it = numerator/denomTarg;
        double ir = numerator/denomRef;
        System.out.println("It " + it);
        System.out.println("Ir " + ir);
        
        System.out.println("Ratio average " + ir/it);
        
        return ir/it;
    }
    private double getDenomTarg(){
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
        if(Math.abs(integrandTargetDenom(xValue, yValue)) >= tol ){
            System.out.println("DenomTarg - increase your ranges xStart, yStart.");
        }
            
        yValue = yEnd;
        if(Math.abs(integrandTargetDenom(xValue, yValue)) >= tol ){
            System.out.println("DenomTarg - increase your ranges xStart, yEnd.");
        }
        xValue = xEnd;
        yValue = yStart;
        if(Math.abs(integrandTargetDenom(xValue, yValue)) >= tol ){
            System.out.println("DenomTarg - increase your ranges xEnd, yStart.");
        }
        
        yValue = yEnd;
        if(Math.abs(integrandTargetDenom(xValue, yValue)) >= tol ){
            System.out.println("DenomTarg - increase your ranges xEnd, yEnd.");
        }
        
        //Now do the "edge" stuff, where either the x or y value, but not both,
        //  are the start or end value
        xValue = xStart;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            if(Math.abs(integrandTargetDenom(xValue, yValue)) >= tol ){
                System.out.println("DenomTarg - increase your ranges xStart, yEdge.");
            }
        }
        xValue = xEnd;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            if(Math.abs(integrandTargetDenom(xValue, yValue)) >= tol ){
                System.out.println("DenomTarg - increase your ranges xEnd, yEdge.");
            }
        }
        yValue = yStart;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(integrandTargetDenom(xValue, yValue)) >= tol ){
                System.out.println("DenomTarg - increase your ranges xEdge, yStart.");
            }
        }
        yValue = yEnd;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(integrandTargetDenom(xValue, yValue)) >= tol ){
                System.out.println("DenomTarg - increase your ranges xEdge, yEnd.");
            }
        }
        
        
        //The actual calculation of the integral.
        //  We have already eliminated the edges' and endpoints' contributions,
        //  because they have been confirmed to be zero by the above code.
        for(int i = 1; i < xN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            for (int j = 1; j < yN; j++) {
                yValue = yStart + j *(yEnd - yStart) / yN;
                    total += 4 * integrandTargetDenom(xValue, yValue);
            }
        }
        //nan could also code the above with xValue = xStart; xValue  += (xEnd - xStart) / xN;
        
        //Now we do the prefix thing
        double prefix = (xEnd - xStart)*(yEnd - yStart) / (4 * xN * yN);
        total *= prefix;
        System.out.println("DenomTarg = " + total);
        return total;
    }
    private double getDenomRef(){
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
        if(Math.abs(integrandRefDenom(xValue, yValue)) >= tol ){
            System.out.println("DenomRef - increase your ranges xStart, yStart.");
        }
            
        yValue = yEnd;
        if(Math.abs(integrandRefDenom(xValue, yValue)) >= tol ){
            System.out.println("DenomRef - increase your ranges xStart, yEnd.");
        }
        xValue = xEnd;
        yValue = yStart;
        if(Math.abs(integrandRefDenom(xValue, yValue)) >= tol ){
            System.out.println("DenomRef - increase your ranges xEnd, yStart.");
        }
        
        yValue = yEnd;
        if(Math.abs(integrandRefDenom(xValue, yValue)) >= tol ){
            System.out.println("DenomRef - increase your ranges xEnd, yEnd.");
        }
        
        //Now do the "edge" stuff, where either the x or y value, but not both,
        //  are the start or end value
        xValue = xStart;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            if(Math.abs(integrandRefDenom(xValue, yValue)) >= tol ){
                System.out.println("DenomRef - increase your ranges xStart, yEdge.");
            }
        }
        xValue = xEnd;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            if(Math.abs(integrandRefDenom(xValue, yValue)) >= tol ){
                System.out.println("DenomRef - increase your ranges xEnd, yEdge.");
            }
        }
        yValue = yStart;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(integrandRefDenom(xValue, yValue)) >= tol ){
                System.out.println("DenomRef - increase your ranges xEdge, yStart.");
            }
        }
        yValue = yEnd;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(integrandRefDenom(xValue, yValue)) >= tol ){
                System.out.println("DenomRef - increase your ranges xEdge, yEnd.");
            }
        }
        
        
        //The actual calculation of the integral.
        //  We have already eliminated the edges' and endpoints' contributions,
        //  because they have been confirmed to be zero by the above code.
        for(int i = 1; i < xN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            for (int j = 1; j < yN; j++) {
                yValue = yStart + j *(yEnd - yStart) / yN;
                    total += 4 * integrandRefDenom(xValue, yValue);
            }
        }
        //nan could also code the above with xValue = xStart; xValue  += (xEnd - xStart) / xN;
        
        //Now we do the prefix thing
        double prefix = (xEnd - xStart)*(yEnd - yStart) / (4 * xN * yN);
        total *= prefix;
        System.out.println("DenomRef = " + total);
        return total;
    }
    private double getNumerator(){
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
        if(Math.abs(integrandNumerator(xValue, yValue)) >= tol ){
            System.out.println("Numerator - increase your ranges xStart, yStart.");
        }
            
        yValue = yEnd;
        if(Math.abs(integrandNumerator(xValue, yValue)) >= tol ){
            System.out.println("Numerator - increase your ranges xStart, yEnd.");
        }
        xValue = xEnd;
        yValue = yStart;
        if(Math.abs(integrandNumerator(xValue, yValue)) >= tol ){
            System.out.println("Numerator - increase your ranges xEnd, yStart.");
        }
        
        yValue = yEnd;
        if(Math.abs(integrandNumerator(xValue, yValue)) >= tol ){
            System.out.println("Numerator - increase your ranges xEnd, yEnd.");
        }
        
        //Now do the "edge" stuff, where either the x or y value, but not both,
        //  are the start or end value
        xValue = xStart;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            if(Math.abs(integrandNumerator(xValue, yValue)) >= tol ){
                System.out.println("Numerator - increase your ranges xStart, yEdge.");
            }
        }
        xValue = xEnd;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            if(Math.abs(integrandNumerator(xValue, yValue)) >= tol ){
                System.out.println("Numerator - increase your ranges xEnd, yEdge.");
            }
        }
        yValue = yStart;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(integrandNumerator(xValue, yValue)) >= tol ){
                        System.out.println("Numerator - increase your ranges xEdge, yStart.");
            }
        }
        yValue = yEnd;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            if(Math.abs(integrandNumerator(xValue, yValue)) >= tol ){
                System.out.println("Numerator - increase your ranges xEdge, yEnd.");
            }
        }
        
        
        //The actual calculation of the integral.
        //  We have already eliminated the edges' and endpoints' contributions,
        //  because they have been confirmed to be zero by the above code.
        for(int i = 1; i < xN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            for (int j = 1; j < yN; j++) {
                yValue = yStart + j *(yEnd - yStart) / yN;
                    total += 4 * integrandNumerator(xValue, yValue);
            }
        }
        //nan could also code the above with xValue = xStart; xValue  += (xEnd - xStart) / xN;
        
        //Now we do the prefix thing
        double prefix = (xEnd - xStart)*(yEnd - yStart) / (4 * xN * yN);
        total *= prefix;
        System.out.println("Numerator = " + total);
        return total;
    }
    private boolean overlap2(double etaReal){
        boolean isOverlap = false;
        
        //This pile o' stuff here calculates the x position of each atom.
        for(int i = 0; i < nAtomsR; i++){
            uR[i] = 0.0;
        }
        
        for(int iCell = 0; iCell < rCells.length; iCell++ ){
            BasisCell cell = rCells[iCell];
            
            double kR = wvR[1].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            //NB all eigenvectors are one, so the term is dropped.
            uR[iCell] += sqrtWvcR[1] * (etaReal * coskR);
            uR[iCell] *= normalizationR;
            atomLocsR[iCell] = uR[iCell] + x0PosR[iCell];
            
        }//end of iCell loop
        
        for (int i=0; i < atomLocsR.length-1; i++){
            if ( (atomLocsR[i] + 0.5) >= (atomLocsR[i+1] - 0.5) ){
                isOverlap = true;
                break;
            }
        }
        
        double repeat = boxR.getBoundary().getBoxSize().getX(0);
        if( (atomLocsR[0]+repeat-0.5) <= (atomLocsR[atomLocsR.length-1]+0.5) ){
            isOverlap = true;
        }
        
        return isOverlap;
    }
    private boolean overlap3(double etaReal, double etaImag){
        boolean isOverlap = false;
        
        //This pile o' stuff here calculates the x position of each atom.
        for(int i = 0; i < nAtomsT; i++){
            uT[i] = 0.0;
        }
        
        for(int iCell = 0; iCell < tCells.length; iCell++ ){
            BasisCell cell = tCells[iCell];
            
            double kR = wvT[1].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);

            //NB all eigenvectors are one, so the term is dropped.
            uT[iCell] += sqrtWvcT[1] * (etaReal * coskR - etaImag * sinkR);
            
            uT[iCell] *= normalizationT;
            atomLocsT[iCell] = uT[iCell] + x0PosT[iCell];
            
        }//end of iCell loop
        
        for (int i=0; i < atomLocsT.length-1; i++){
            if ( (atomLocsT[i] + 0.5) >= (atomLocsT[i+1] - 0.5) ){
                isOverlap = true;
                break;
            }
        }
        
        double repeat = boxT.getBoundary().getBoxSize().getX(0);
        if( (atomLocsT[0]+repeat-0.5) <= (atomLocsT[atomLocsT.length-1]+0.5) ){
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
    private double integrandRefDenom(double x, double y){
        double value;
        
        //Hard Rod + harmonic, 2 rods
        if (overlap2(x)) {
            value = 0.0;
        }else{
            value = Math.exp(-0.5 * omega2R[1][0] * y*y);
        }
        return value;
    }
    private double integrandTargetDenom(double x, double y){
        double value;
        
//      //Hard rod N = 3
        if (overlap3(x, y)) {
            value = 0;
        }else{
            value = 1;
        }
        return value;
    }
    private double integrandNumerator(double x, double y){
        double value;
        
        if (!overlap2(x) && !overlap3(x, y)){
            double dork = Math.exp(-0.5 * omega2R[1][0] * y*y);
            value = dork / (1 + alpha * dork);
        } else {
            value = 0;
        }
        return value;
    }
    
    public static void main(String[] args) {
      double xStart = -1.0;
      double yStart = -1.0;
      double xEnd = 1.0;
      double yEnd = 1.0;
      
      int xN = 400;
      int yN = 400;
      
      int nAtomsTarget = 3;
      double density = 0.7;
      double alpha = 1.5;
      
      System.out.println("Alpha " + alpha);
      DoubleIntegralEtas di = new DoubleIntegralEtas(nAtomsTarget, density, alpha);
      di.setIntegrationParameters(xStart, xEnd, yStart, yEnd, xN, yN);
      di.calculate();
    }
}
