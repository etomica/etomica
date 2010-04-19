package etomica.models.oneDHardRods;

import etomica.api.IAtomType;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
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
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Null;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume a mode & rod 
 * are removed
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImageSubtract extends DataSourceScalar {

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
    private double[][] simOmegaSquared, oneOverOmega2;

    protected final IRandom random;
    private IBox box;
    private int numAtoms;
    private Boundary bdry;
    private NormalModes nm;
    WaveVectorFactory waveVectorFactory;
    private double etas[];
    private int min, max;
    
    public MeterDifferentImageSubtract(ISimulation sim, ISpace space, double temp, 
            CoordinateDefinition simCD, NormalModes simNM, IBox otherBox){
        super("MeterSubtract", Null.DIMENSION);
        this.random = sim.getRandom();
        this.temperature = temp;
        
        simWaveVectors = simNM.getWaveVectorFactory().getWaveVectors();
        this.simCDef = simCD;
        simCDim = simCD.getCoordinateDim();
        simEigenVectors = simNM.getEigenvectors();
        simWVCoeff = simNM.getWaveVectorFactory().getCoefficients();
        simRealT = new double[simCDim];
        simImagT = new double[simCDim];
        double[][] tempO2 = simNM.getOmegaSquared();
        for(int i = 0; i < tempO2.length; i++ ){
            for(int j = 0; j < tempO2[0].length; j++){
                simOmegaSquared[i][j] = Math.sqrt(tempO2[i][j]);
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
        
        //nan This will change when we have more than one atoms per cell
        int[] nCells = new int[]{numAtoms};
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
        tempO2 = nm.getOmegaSquared();
        for(int i = 0; i < tempO2.length; i++ ){
            for(int j = 0; j < tempO2[0].length; j++){
                oneOverOmega2[i][j] = 1/Math.sqrt(tempO2[i][j]);
            }
        }
        wvCoeff = nm.getWaveVectorFactory().getCoefficients();
        sqrtWVC = new double[wvCoeff.length];
        for (int i =0; i < wvCoeff.length; i++){
            sqrtWVC[i] = Math.sqrt(2*wvCoeff[i]);
        }
        
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
        
        etas = new double[2 * simWaveVectors.length * simCDim];
        min = space.D() * (numAtoms - 1); 
        max = space.D() * (simCDef.getBox().getLeafList().getAtomCount() - 1) ;
    }
    
    public double getDataAsScalar() {
//        IAtomList atomlist = box.getLeafList();
//        for (int i = 0; i < atomlist.getAtomCount(); i++){
//            System.out.println("start i " + atomlist.getAtom(i).getPosition().getX(0));
//        }
        BasisCell[] cells = cDef.getBasisCells();
        BasisCell cell = cells[0];
        //nan this makes it 1D
        newU = new double[cDim];
        for(int i=0; i < newU.length; i++){
            newU[i] = 0.0;
        }
        
        //Calculate normal mode coordinates of simulation system.
        double[] realCoord = new double[simWaveVectors.length];
        double[] imagCoord = new double[simWaveVectors.length];
        int etaCount = 0;
        for (int iWV = 0; iWV < simWaveVectors.length; iWV++){
            simCDef.calcT(simWaveVectors[iWV], simRealT, simImagT);
            realCoord[iWV] = 0.0;
            imagCoord[iWV] = 0.0;
            for (int i = 0; i < simCDim; i++){
                for (int j = 0; j < simCDim; j++){
                    realCoord[iWV] += simEigenVectors[iWV][i][j] * simRealT[j];
                    imagCoord[iWV] += simEigenVectors[iWV][i][j] * simImagT[j];
                }
            }
            if(simWVCoeff[iWV] == 1.0){
                realCoord[iWV] *= Math.sqrt(2);
                imagCoord[iWV] *= Math.sqrt(2);
                //nan omega2[wv][evect]
                etas[etaCount] = realCoord[iWV] * simOmegaSquared[iWV][0];
                etaCount++;
                etas[etaCount] = imagCoord[iWV] * simOmegaSquared[iWV][0];
                etaCount++;
            } else {
                etas[etaCount] = realCoord[iWV] * simOmegaSquared[iWV][0];
                etaCount++;
            }
        }
        
        //CALCULATE A HARMONIC ENERGY FOR THE SUBTRACTED NORMAL MODES.
        //nan THERE IS NO ORGANIZATION TO WHICH ETAS ARE ZEROED!!!!!!
        
        //Calculate harmonic, and zero out the subtracted normal modes.
        double harmonic = 0.0;
        for (int i = min; i < max; i++){
            //The omega^2 term is dropped because of the scaling.
            harmonic = 0.5 * etas[i] * etas[i];
            etas[i] = 0.0;
        }
        
        //Calculate the positions for the meter's system
        for (int iCell = 0; iCell < cells.length; iCell++){
            cell = cells[iCell];
            for (int j = 0; j < cDim; j++) {
                newU[j] = 0.0;
            }
            etaCount = 0;
            for (int iWV = 0; iWV < waveVectors.length; iWV++){
                //Calculate the change in positions.
                double kR = waveVectors[iWV].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int iMode = 0; iMode < cDim; iMode++){
                    for (int iCD = 0; iCD < cDim; iCD++){
                        if (wvCoeff[iWV] == 1.0){
                            newU[iCD] += sqrtWVC[iWV] 
                                * eigenVectors[iWV][iMode][iCD] * oneOverOmega2[iWV][0]
                                * (etas[etaCount] * coskR - etas[etaCount+1] * sinkR);
                            etaCount += 2;
                        } else {
                            newU[iCD] += sqrtWVC[iWV] 
                                * eigenVectors[iWV][iMode][iCD] * oneOverOmega2[iWV][0]
                                * (etas[etaCount] * coskR);
                            etaCount ++;
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
        
//        return Math.exp(-1*(meterPE.getDataAsScalar()+harmonic));
        return meterPE.getDataAsScalar() + harmonic;
    }

    
}
