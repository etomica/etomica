package etomica.models.oneDHardRods;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.NormalModes;
import etomica.normalmode.WaveVectorFactory;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
import etomica.units.Null;
import etomica.normalmode.BasisBigCell;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume a mode & rod 
 * are removed
 * 
 * Should only be used when one system is exactly twice the size of the other
 * system.
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImageSubtractDoubleSize extends DataSourceScalar {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MeterDifferentImageSubtract";
    
    public int nInsert, counter;
    private MeterPotentialEnergy meterPE;
    private CoordinateDefinition cDef, simCDef;
    private int cDim, simCDim;
    private IVectorMutable[] waveVectors, simWaveVectors;
    private double[] simRealT, simImagT;
    protected double temperature;
    private double[] newU;
    private double[] wvCoeff, simWVCoeff, sqrtWVC;
    private double[][][] eigenVectors, simEigenVectors;
    private double[][] sqrtSimOmega2, oneOverSqrtOmega2;

    protected final IRandom random;
    private IBox box;
    private int numAtoms;
    private IBoundary bdry;
    private NormalModes nm;
    WaveVectorFactory waveVectorFactory;
    private double etas[];
    private int maxEta;
    private double scaling;
    
    
    public MeterDifferentImageSubtractDoubleSize(ISimulation sim, ISpace space, 
            CoordinateDefinition simCD, NormalModes simNM, CoordinateDefinition
            otherCD, PotentialMasterList potentialMaster, int[] otherNCells, 
            NormalModes otherNM){
        this(sim, space, simCD, simNM, otherCD, potentialMaster, 
                otherNCells, otherNM, "file");
    }
    public MeterDifferentImageSubtractDoubleSize(ISimulation sim, ISpace space, 
            CoordinateDefinition simCD, NormalModes simNM, CoordinateDefinition
            otherCD, PotentialMasterList potentialMaster, int[] otherNCells, 
            NormalModes otherNM, String otherFilename){
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
        for(int i = 0; i < tempO2.length; i++ ){
            for(int j = 0; j < tempO2[0].length; j++){
                if(Double.isInfinite(tempO2[i][j])){
                    sqrtSimOmega2[i][j] = tempO2[i][j];
                } else {
                    sqrtSimOmega2[i][j] = Math.sqrt(tempO2[i][j]);
                }
            }
        }
        
        double density = simCDef.getBox().getLeafList().getAtomCount() / 
            simCDef.getBox().getBoundary().volume();
        numAtoms = otherCD.getBox().getLeafList().getAtomCount();
        box = new Box(space);
        sim.addBox(box);
        box.setNMolecules(sim.getSpecies(0), numAtoms);
        
        if (space.D() == 1) {
            bdry = new BoundaryRectangularPeriodic(space, numAtoms/ density);
        } else {
            bdry = new BoundaryRectangularPeriodic(space, 1.0);
            IVector edges = otherCD.getBox().getBoundary().getBoxSize();
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
        for (int i =0; i < wvCoeff.length; i++){
            sqrtWVC[i] = Math.sqrt(2*wvCoeff[i]);
        }
        tempO2 = nm.getOmegaSquared();
        oneOverSqrtOmega2 = new double[tempO2.length][tempO2[0].length];
        for(int i = 0; i < tempO2.length; i++ ){
            for(int j = 0; j < tempO2[0].length; j++){
                oneOverSqrtOmega2[i][j] = 1/Math.sqrt(tempO2[i][j]);
            }
        }
        
        //Create the scaling factor
        scaling = 0.0;
        for(int iWV = 0; iWV < simWaveVectors.length; iWV++){
            for(int iMode = 0; iMode < simCDim; iMode++){
                if(!Double.isInfinite(sqrtSimOmega2[iWV][iMode])){
                    scaling += Math.log(sqrtSimOmega2[iWV][iMode]);
                    if(simWVCoeff[iWV] == 1.0){
                        scaling += Math.log(sqrtSimOmega2[iWV][iMode]);
                    }
                }
            }
        }
        for (int iWV = 0; iWV < waveVectors.length; iWV++){
            for(int iMode = 0; iMode < cDim; iMode++){
                if(!(oneOverSqrtOmega2[iWV][iMode] == 0.0)){
                    scaling += Math.log(oneOverSqrtOmega2[iWV][iMode]);
                    if (wvCoeff[iWV] == 1.0){
                        scaling += Math.log(oneOverSqrtOmega2[iWV][iMode]);
                    }
                }
            }
        }
        
        potentialMaster.getNeighborManager(box).reset();
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        
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
        double[] realCoord = new double[simCDim];
        double[] imagCoord = new double[simCDim];
        
        simCDef.calcT(simWaveVectors[1], simRealT, simImagT);
        for (int iMode = 0; iMode < simCDim; iMode++){
            realCoord[iMode] = 0.0;
            imagCoord[iMode] = 0.0;
            for (int j = 0; j < simCDim; j++){
                realCoord[iMode] += simEigenVectors[1][iMode][j] * simRealT[j];
                imagCoord[iMode] += simEigenVectors[1][iMode][j] * simImagT[j];
            }
            if(simWVCoeff[1] == 1.0){
                realCoord[iMode] *= Math.sqrt(2);
                imagCoord[iMode] *= Math.sqrt(2);
            }
        }

        
        //CALCULATE A HARMONIC ENERGY FOR THE SUBTRACTED NORMAL MODES.
        
        //Calculate harmonic, and zero out the subtracted normal modes.
        double harmonic = 0.0;
        if(simWVCoeff[1] == 1.0){
            for (int iMode = 0; iMode < simCDim; iMode++){
                if(!(sqrtSimOmega2[1][iMode] == Double.POSITIVE_INFINITY)){
                    harmonic += 0.5 * realCoord[iMode] * sqrtSimOmega2[1][iMode]
                        * realCoord[iMode] * sqrtSimOmega2[1][iMode];
                    harmonic += 0.5 * imagCoord[iMode] * sqrtSimOmega2[1][iMode]
                        * imagCoord[iMode] * sqrtSimOmega2[1][iMode];
                }
            }
        } else {
            for (int iMode = 0; iMode < simCDim; iMode++){
                if(!(sqrtSimOmega2[1][iMode] == Double.POSITIVE_INFINITY)){
                    harmonic += 0.5 * realCoord[iMode] * sqrtSimOmega2[1][iMode]
                        * realCoord[iMode] * sqrtSimOmega2[1][iMode];
                }
            }
        }
        
        //Calculate the positions for the meter's system
//        for (int iCell = 0; iCell < cells.length; iCell++){
//            cell = cells[iCell];
//            for (int j = 0; j < cDim; j++) {
//                newU[j] = 0.0;
//            }
//            etaCount = 0;   //etaCount counts through "wv" for etas.
//            for (int iWV = 0; iWV < waveVectors.length; iWV++){
//                //Calculate the change in positions.
//                double kR = waveVectors[iWV].dot(cell.cellPosition);
//                double coskR = Math.cos(kR);
//                double sinkR = Math.sin(kR);
//                if(wvCoeff[iWV] == 1.0){
//                    for (int iMode = 0; iMode < cDim; iMode++){
//                        if(etaCount == maxEta) {break;}
//                        if(!(oneOverSqrtOmega2[iWV][iMode] == 0.0)){
//                            for (int iCD = 0; iCD < cDim; iCD++){
//                                newU[iCD] += sqrtWVC[iWV] * eigenVectors[iWV][iMode][iCD] 
//                                    * oneOverSqrtOmega2[iWV][iMode]
//                                    * (etas[etaCount] * coskR - etas[etaCount+1]* sinkR);
//                            }
//                            etaCount += 2;
//                        }
//                    }
//                } else {
//                    for (int iMode = 0; iMode < cDim; iMode++){
//                        if(etaCount == maxEta) {break;}
//                        if(!(oneOverSqrtOmega2[iWV][iMode] == 0.0)){
//                            for (int iCD = 0; iCD < cDim; iCD++){
//                                newU[iCD] += sqrtWVC[iWV] * eigenVectors[iWV][iMode][iCD] 
//                                    * oneOverSqrtOmega2[iWV][iMode]
//                                    * (etas[etaCount] * coskR);
//                            }
//                            etaCount += 1;
//                        }
//                    }
//                }
//            }
//            
//            if(etaCount != maxEta){
//                throw new IllegalStateException("etaCount != maxEta in MeterDifferentImageSubtract");
//            }
//            
//            for (int i=0; i<cDim; i++) {
//                newU[i] *= normalization;
//            }
//            cDef.setToU(cells[iCell].molecules, newU);
//        }
        
        return meterPE.getDataAsScalar() + harmonic;
    }


    public IBox getBox(){
        return box;
    }
    
    /**
     * @return the natural logarithm of the scaling
     */
    public double getScaling() {
        return scaling;
    }
}
