/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.space.Boundary;
import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.NormalModes;
import etomica.normalmode.WaveVectorFactory;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.units.dimensions.Null;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume a mode & rod 
 * are removed
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImageSubtract extends DataSourceScalar {
    protected static final long serialVersionUID = 1L;
    protected static final String APP_NAME = "MeterDifferentImageSubtract";
    
    public int nInsert, counter;
    protected MeterPotentialEnergy meterPE;
    protected CoordinateDefinition cDef, simCDef;
    protected int cDim, simCDim;
    protected Vector[] waveVectors, simWaveVectors;
    protected double[] simRealT, simImagT;
    protected double temperature;
    protected double[] newU;
    protected double[] wvCoeff, simWVCoeff, sqrtWVC;
    protected double[][][] eigenVectors, simEigenVectors;
    protected double[][] sqrtSimOmega2, oneOverSqrtOmega2;

    protected final IRandom random;
    protected Box box;
    protected int numAtoms;
    protected Boundary bdry;
    protected NormalModes nm;
    WaveVectorFactory waveVectorFactory;
    protected double etas[];
    protected int maxEta;
    protected double scaling;
    
    
    public MeterDifferentImageSubtract(Simulation sim, Space space,
                                       CoordinateDefinition simCD, NormalModes simNM, CoordinateDefinition
            otherCD, PotentialMasterList potentialMaster, int[] otherNCells,
                                       NormalModes otherNM){
        this(sim, space, simCD, simNM, otherCD, potentialMaster, 
                otherNCells, otherNM, "file");
    }
    public MeterDifferentImageSubtract(Simulation sim, Space space,
                                       CoordinateDefinition simCD, NormalModes simNM, CoordinateDefinition
            otherCD, PotentialMasterList potentialMaster, int[] otherNCells,
                                       NormalModes otherNM, String otherFilename) {
        super("MeterSubtract", Null.DIMENSION);
        this.random = sim.getRandom();

        simWaveVectors = simNM.getWaveVectorFactory().getWaveVectors();
        this.simCDef = simCD;
        simCDim = simCD.getCoordinateDim();
        cDim = otherCD.getCoordinateDim();
        simEigenVectors = simNM.getEigenvectors();
        simWVCoeff = simNM.getWaveVectorFactory().getCoefficients();
        simRealT = new double[simCDim];
        simImagT = new double[simCDim];
        double[][] tempO2 = simNM.getOmegaSquared();
        sqrtSimOmega2 = new double[tempO2.length][tempO2[0].length];
        for (int i = 0; i < tempO2.length; i++) {
            for (int j = 0; j < tempO2[0].length; j++) {
                if (Double.isInfinite(tempO2[i][j])) {
                    sqrtSimOmega2[i][j] = tempO2[i][j];
                } else {
                    sqrtSimOmega2[i][j] = Math.sqrt(tempO2[i][j]);
                }
            }
        }

        double density = simCDef.getBox().getLeafList().getAtomCount() /
                simCDef.getBox().getBoundary().volume();
        numAtoms = otherCD.getBox().getLeafList().getAtomCount();
        box = sim.makeBox();
        box.setNMolecules(sim.getSpecies(0), numAtoms);

        if (space.D() == 1) {
            bdry = new BoundaryRectangularPeriodic(space, numAtoms / density);
        } else {
            bdry = new BoundaryRectangularPeriodic(space, 1.0);
            Vector edges = otherCD.getBox().getBoundary().getBoxSize();
            bdry.setBoxSize(edges);
        }
        box.setBoundary(bdry);

        cDef = new CoordinateDefinitionLeaf(box, otherCD.getPrimitive(),
                otherCD.getBasis(), space);
        cDef.initializeCoordinates(otherNCells);


        nm = otherNM;
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        waveVectors = nm.getWaveVectorFactory().getWaveVectors();
        eigenVectors = nm.getEigenvectors();

        wvCoeff = nm.getWaveVectorFactory().getCoefficients();
        sqrtWVC = new double[wvCoeff.length];
        for (int i = 0; i < wvCoeff.length; i++) {
            sqrtWVC[i] = Math.sqrt(2 * wvCoeff[i]);
        }
        tempO2 = nm.getOmegaSquared();
        oneOverSqrtOmega2 = new double[tempO2.length][tempO2[0].length];
        for (int i = 0; i < tempO2.length; i++) {
            for (int j = 0; j < tempO2[0].length; j++) {
                oneOverSqrtOmega2[i][j] = 1 / Math.sqrt(tempO2[i][j]);
            }
        }

        //Create the scaling factor
        scaling = 0.0;
        for (int iWV = 0; iWV < simWaveVectors.length; iWV++) {
            for (int iMode = 0; iMode < simCDim; iMode++) {
                if (!Double.isInfinite(sqrtSimOmega2[iWV][iMode])) {
                    scaling += Math.log(sqrtSimOmega2[iWV][iMode]);
                    if (simWVCoeff[iWV] == 1.0) {
                        scaling += Math.log(sqrtSimOmega2[iWV][iMode]);
                    }
                }
            }
        }
        for (int iWV = 0; iWV < waveVectors.length; iWV++) {
            for (int iMode = 0; iMode < cDim; iMode++) {
                if (!(oneOverSqrtOmega2[iWV][iMode] == 0.0)) {
                    scaling += Math.log(oneOverSqrtOmega2[iWV][iMode]);
                    if (wvCoeff[iWV] == 1.0) {
                        scaling += Math.log(oneOverSqrtOmega2[iWV][iMode]);
                    }
                }
            }
        }

        potentialMaster.getNeighborManager(box).reset();
        meterPE = new MeterPotentialEnergy(potentialMaster, box);

        etas = new double[space.D() * (simCDef.getBox().getLeafList().getAtomCount() - 1)];
        maxEta = space.D() * (numAtoms - 1);

    }
    
    public double getDataAsScalar() {
        BasisCell[] simCells = simCDef.getBasisCells();
        BasisCell[] cells = cDef.getBasisCells();
        double normalization = 1/Math.sqrt(cells.length);
        BasisCell cell = simCells[0];
        newU = new double[cDim];
        for(int i=0; i < newU.length; i++){
            newU[i] = 0.0;
        }
        
        //Calculate normal mode coordinates of simulation system.
        // First index is wavevector, second index is mode.
        double[][] realCoord = new double[simWaveVectors.length][simCDim];
        double[][] imagCoord = new double[simWaveVectors.length][simCDim];
        for (int iWV = 0; iWV < simWaveVectors.length; iWV++){
            simCDef.calcT(simWaveVectors[iWV], simRealT, simImagT);
            for (int iMode = 0; iMode < simCDim; iMode++){
                realCoord[iWV][iMode] = 0.0;
                imagCoord[iWV][iMode] = 0.0;
                for (int j = 0; j < simCDim; j++){
                    realCoord[iWV][iMode] += simEigenVectors[iWV][iMode][j] * simRealT[j];
                    imagCoord[iWV][iMode] += simEigenVectors[iWV][iMode][j] * simImagT[j];
                }
                if(simWVCoeff[iWV] == 1.0){
                    realCoord[iWV][iMode] *= Math.sqrt(2);
                    imagCoord[iWV][iMode] *= Math.sqrt(2);
                }
            }
        }
        
        //Scale and transfer the normal mode coordinates to etas.
        int etaCount = 0;
        for (int iWV = 0; iWV < simWaveVectors.length; iWV++){
            if(simWVCoeff[iWV] == 1.0){
                for (int iMode = 0; iMode < simCDim; iMode++){
                    if(!(sqrtSimOmega2[iWV][iMode] == Double.POSITIVE_INFINITY)){
                        etas[etaCount] = realCoord[iWV][iMode] * 
                                sqrtSimOmega2[iWV][iMode];
                        etas[etaCount+1] = imagCoord[iWV][iMode] * 
                                sqrtSimOmega2[iWV][iMode];
                        etaCount +=2;
                    }
                }
            } else {
                for (int iMode = 0; iMode < simCDim; iMode++){
                    if(!(sqrtSimOmega2[iWV][iMode] == Double.POSITIVE_INFINITY)){
                        etas[etaCount] = realCoord[iWV][iMode] * 
                                sqrtSimOmega2[iWV][iMode];
                        etaCount++;
                    }
                }
            }
        }
        
        //CALCULATE A HARMONIC ENERGY FOR THE SUBTRACTED NORMAL MODES.
        
        //Calculate harmonic, and zero out the subtracted normal modes.
        double harmonic = 0.0;
        for (int i = maxEta; i < etas.length; i++){
            //The omega^2 term is dropped because of the scaling.
            harmonic += 0.5 * etas[i] * etas[i];
        }
        
        //Calculate the positions for the meter's system
        for (int iCell = 0; iCell < cells.length; iCell++){
            cell = cells[iCell];
            for (int j = 0; j < cDim; j++) {
                newU[j] = 0.0;
            }
            etaCount = 0;   //etaCount counts through "wv" for etas.
            for (int iWV = 0; iWV < waveVectors.length; iWV++){
                //Calculate the change in positions.
                double kR = waveVectors[iWV].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                if(wvCoeff[iWV] == 1.0){
                    for (int iMode = 0; iMode < cDim; iMode++){
                        if(etaCount == maxEta) {break;}
                        if(!(oneOverSqrtOmega2[iWV][iMode] == 0.0)){
                            for (int iCD = 0; iCD < cDim; iCD++){
                                newU[iCD] += sqrtWVC[iWV] * eigenVectors[iWV][iMode][iCD] 
                                    * oneOverSqrtOmega2[iWV][iMode]
                                    * (etas[etaCount] * coskR - etas[etaCount+1]* sinkR);
                            }
                            etaCount += 2;
                        }
                    }
                } else {
                    for (int iMode = 0; iMode < cDim; iMode++){
                        if(etaCount == maxEta) {break;}
                        if(!(oneOverSqrtOmega2[iWV][iMode] == 0.0)){
                            for (int iCD = 0; iCD < cDim; iCD++){
                                newU[iCD] += sqrtWVC[iWV] * eigenVectors[iWV][iMode][iCD] 
                                    * oneOverSqrtOmega2[iWV][iMode]
                                    * (etas[etaCount] * coskR);
                            }
                            etaCount += 1;
                        }
                    }
                }
            }
            
            if(etaCount != maxEta){
                throw new IllegalStateException("etaCount != maxEta in MeterDifferentImageSubtract");
            }
            
            for (int i=0; i<cDim; i++) {
                newU[i] *= normalization;
            }
            cDef.setToU(cells[iCell].molecules, newU);
        }
        
        
        return meterPE.getDataAsScalar() + harmonic;
    }


    public Box getBox(){
        return box;
    }
    
    /**
     * @return the natural logarithm of the scaling
     */
    public double getScaling() {
        return scaling;
    }
}
