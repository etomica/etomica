package etomica.models.oneDHardRods;

import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.NormalModes;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Null;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume an extra mode & rod 
 * are added
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImageAdd extends DataSourceScalar {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MeterDifferentImageAdd";
    
    public int nInsert, counter;
    private MeterPotentialEnergy meterPE;
    private CoordinateDefinition cDef, simCDef;
    private int cDim, spaceD;
    private IVectorMutable[] waveVectors, simWaveVectors;
    private double[] simRealT, simImagT;
    protected double temperature;
    private double[] newU;
    private double[] wvCoeff, simWVCoeff, sqrtWVC;
    private double[][] oneOverSqrtOmega2, sqrtSimOmega2; //These are already made sqrt.
    private double[][][] eigenVectors, simEigenVectors;
    double[] gaussCoord;
    
    protected final IRandom random;
    public IBox box;
    private int numAtoms;
    private IBoundary bdry;
    private NormalModes nm;
    WaveVectorFactory waveVectorFactory;
    private double etas[];
    private double scaling;
    private int numberOfReal, numberOfImag;
    
    
    public MeterDifferentImageAdd(ISimulation sim, ISpace space, double temp, 
            CoordinateDefinition simCD, NormalModes simNM, IBox otherBox){
        
        super("MeterAdd", Null.DIMENSION);
        this.random = sim.getRandom();
        this.temperature = temp;
        
        simWaveVectors = simNM.getWaveVectorFactory().getWaveVectors();
        this.simCDef = simCD;
        cDim = simCD.getCoordinateDim();
        simEigenVectors = simNM.getEigenvectors();
        simWVCoeff = simNM.getWaveVectorFactory().getCoefficients();
        simRealT = new double[cDim];
        simImagT = new double[cDim];
        double[][] omegaTemp = simNM.getOmegaSquared();
        sqrtSimOmega2 = new double[omegaTemp.length][omegaTemp[0].length];
        for(int i = 0; i < omegaTemp.length; i++){
            for(int j = 0; j < omegaTemp[0].length; j++){
                if(Double.isInfinite(omegaTemp[i][j])){
                    sqrtSimOmega2[i][j] = omegaTemp[i][j];
                } else {
                    sqrtSimOmega2[i][j] = Math.sqrt(omegaTemp[i][j]);
                }
            }
        }
        
        double density = simCDef.getBox().getLeafList().getAtomCount() / 
                simCDef.getBox().getBoundary().volume();
        numAtoms = otherBox.getLeafList().getAtomCount();
        box = new Box(space);
        sim.addBox(box);
        box.setNMolecules(sim.getSpecies(0), numAtoms);
        bdry = new BoundaryRectangularPeriodic(space, numAtoms/density);
        box.setBoundary(bdry);
        
        //scale for the number of cells
        int[] nCells = new int[]{numAtoms * simCDef.getBasisCells().length
                / simCDef.getBox().getLeafList().getAtomCount()};
        cDef = new CoordinateDefinitionLeaf(box, simCDef.getPrimitive(), 
                simCDef.getBasis(), space);
        cDef.initializeCoordinates(nCells);
        cDim = cDef.getCoordinateDim();

        nm = new NormalModes1DHR(box.getBoundary(), numAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(temperature);
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        waveVectors = nm.getWaveVectorFactory().getWaveVectors();
        eigenVectors = nm.getEigenvectors();
        wvCoeff = nm.getWaveVectorFactory().getCoefficients();
        sqrtWVC = new double[wvCoeff.length];
        for (int i =0; i < wvCoeff.length; i++){
            sqrtWVC[i] = Math.sqrt(2*wvCoeff[i]);
        }
        omegaTemp = nm.getOmegaSquared();
        scaling = 1.0;
        oneOverSqrtOmega2 = new double[omegaTemp.length][omegaTemp[0].length];
        for (int i = 0; i < oneOverSqrtOmega2.length; i++) {
            for (int j = 0; j < oneOverSqrtOmega2[i].length; j++) {
                oneOverSqrtOmega2[i][j] = Math.sqrt(1.0/(omegaTemp[i][j]));
            }
        }
        
        //nan 1D stuff - will need to change
        PotentialMasterList potentialMaster = new PotentialMasterList(sim, space);
        Potential2 potential = new P2HardSphere(space, 1.0, true);
        potential = new P2XOrder(space, (Potential2HardSpherical)potential);
        potential.setBox(box);
        potentialMaster.addPotential(potential, new IAtomType[] {
                ((SpeciesSpheresMono)sim.getSpecies(0)).getLeafType(), 
                ((SpeciesSpheresMono)sim.getSpecies(0)).getLeafType()});
        double neighborRange = 1.01/density;
        potentialMaster.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMaster.getNeighborManager(box).reset();
        
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        
        etas = new double[space.D() * (numAtoms - 1)];
        
        gaussCoord = new double[space.D() *(numAtoms - simCDef.getBox().getLeafList().getAtomCount())];
        
        spaceD = space.D();
        
        //Create the scaling factor
        scaling = 1.0;
        for(int iWV = 0; iWV < simWaveVectors.length; iWV++){
            for(int iMode = 0; iMode < cDim; iMode++){
                if(!Double.isInfinite(sqrtSimOmega2[iWV][iMode])){
                    scaling *= sqrtSimOmega2[iWV][iMode];
                    if(simWVCoeff[iWV] == 1.0){
                        scaling *= sqrtSimOmega2[iWV][iMode];
                    }
                }
            }
        }
        for(int iWV = 0; iWV < waveVectors.length; iWV++){
            for(int iMode = 0; iMode < cDim; iMode++){
                if(!(oneOverSqrtOmega2[iWV][iMode] == 0.0)){
                    scaling *= oneOverSqrtOmega2[iWV][iMode];
                    if (wvCoeff[iWV] == 1.0){
                        scaling *= oneOverSqrtOmega2[iWV][iMode];
                    }
                }
            }
        }
        
        //create the numbers used for loops later on.
        numberOfReal = 0;
        numberOfImag = 0;
        for(int iWV = 0; iWV < simWaveVectors.length; iWV++){
            for(int iMode = 0; iMode < cDim; iMode++){
                if(!Double.isInfinite(sqrtSimOmega2[iWV][iMode])){
                    numberOfReal++;
                    if (simWVCoeff[iWV] == 1.0){
                        numberOfImag++;
                    }
                }
            }
        }
        
   }
    
    public double getDataAsScalar() {
        
        BasisCell[] simCells = simCDef.getBasisCells();
        BasisCell[] cells = cDef.getBasisCells();
        BasisCell cell = simCells[0];
        //nan this makes it 1D
        newU = new double[cDim];
        for(int i=0; i < newU.length; i++){
            newU[i] = 0.0;
        }
        
        //Calculate normal mode coordinates of simulation system.
        double[] realCoord = new double[waveVectors.length * cDim];
        double[] imagCoord = new double[waveVectors.length * cDim];
        for (int iWV = 0; iWV < simWaveVectors.length; iWV++){
            simCDef.calcT(simWaveVectors[iWV], simRealT, simImagT);
            realCoord[iWV] = 0.0;
            imagCoord[iWV] = 0.0;
            for (int iMode = 0; iMode < cDim; iMode++){
                for (int j = 0; j < cDim; j++){
                    realCoord[iWV] += simEigenVectors[iWV][iMode][j] * simRealT[j];
                    imagCoord[iWV] += simEigenVectors[iWV][iMode][j] * simImagT[j];
                }
            }
            if(simWVCoeff[iWV] == 1.0){
                realCoord[iWV] *= Math.sqrt(2);
                imagCoord[iWV] *= Math.sqrt(2);
            }
        }
        
        //Scale and transfer the normal mode coordinates to etas.
        int realCount = 0;
        int imagCount = 1;
        for(int iWV = 0; iWV < simWaveVectors.length; iWV++){
            for(int iMode = 0; iMode < cDim; iMode++){
                if(!Double.isInfinite(sqrtSimOmega2[iWV][iMode])){
                    etas[realCount] = realCoord[iWV] * sqrtSimOmega2[iWV][iMode];
                    realCount += 2;
                    if(simWVCoeff[iWV] == 1.0){
                        etas[imagCount] = imagCoord[iWV] * sqrtSimOmega2[iWV][iMode];
                        imagCount += 2;
                    }
                }
            }
        }
        
        //Create the last normal mode coordinates from the Gaussian distribution 
        for (int count = realCount; count < etas.length; count+=2){
            gaussCoord[count-realCount] = random.nextGaussian();
            etas[count] = gaussCoord[count-realCount];
        }
        for(int count = imagCount; count < etas.length; count+=2){
            gaussCoord[count-imagCount] = random.nextGaussian();
            etas[count] = gaussCoord[count-imagCount];
        }
        
        //Calculate the positions for the meter's system
        for (int iCell = 0; iCell < cells.length; iCell++){
            cell = cells[iCell];
            for (int j = 0; j < cDim; j++) {
                newU[j] = 0.0;
            }
            int etaCount = 0;   //etaCount counts through "wv" for etas.
            for (int wvcount = 0; wvcount < waveVectors.length; wvcount++){
                //Calculate the change in positions.
                double kR = waveVectors[wvcount].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int iMode = 0; iMode < cDim; iMode++){
                    if(!(oneOverSqrtOmega2[wvcount][iMode] == 0.0)){
                        for (int iCD = 0; iCD < cDim; iCD++){
                            if(etaCount+1 < etas.length){
                               newU[iCD] += sqrtWVC[wvcount] * eigenVectors[wvcount][iMode][iCD] 
                                    * oneOverSqrtOmega2[wvcount][iMode]
                                    * (etas[etaCount] * coskR - etas[etaCount+1]* sinkR);
                               etaCount += 2;
                            } else {
                                newU[iCD] += sqrtWVC[wvcount] * eigenVectors[wvcount][iMode][iCD]
                                    * oneOverSqrtOmega2[wvcount][iMode]
                                    * (etas[etaCount] * coskR);
                                etaCount += 1;
                            }
                        }
                    }
                }
            }
            
            double normalization = 1/Math.sqrt(cells.length);
            for (int i=0; i<cDim; i++) {
                newU[i] *= normalization;
            }
            
            cDef.setToU(cells[iCell].molecules, newU);
        }
        
//        return Math.exp(-1*meterPE.getDataAsScalar());
        return meterPE.getDataAsScalar();
    }
    
    public double[] getGaussian(){
        return gaussCoord;
    }

    public double getScaling() {
        return scaling;
    }
    
}
