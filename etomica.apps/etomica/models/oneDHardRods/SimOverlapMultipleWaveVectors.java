package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.DataPump;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.WaveVectorFactory;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

public class SimOverlapMultipleWaveVectors extends Simulation {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimOverlapABC";
    Primitive primitive;
    int[] nCells;
    NormalModes1DHR nm;
    double bennettParam;       //adjustable parameter - Bennett's parameter
    public IntegratorOverlap integratorSim; //integrator for the whole simulation
    public DataSourceVirialOverlap dsvo;
    public BasisMonatomic basis;
    ActivityIntegrate activityIntegrate;
    
    IntegratorMC[] integrators;
    protected int blockSize;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IEtomicaDataSource[] meters;
    public IBox boxTarget, boxRef;
    public Boundary boundaryTarget, boundaryRef;
    MCMoveChangeMode changeMove;
    MCMoveCompareMultipleModes compareMove;
    MeterPotentialEnergy meterAinB, meterAinA;
    MeterCompareMultipleModesBrute meterBinA, meterBinB;
    MeterCompareTest meterTestBinA, meterTestBinB;
    MeterOverlap meterOverlapInA, meterOverlapInB;
    
    
    public SimOverlapMultipleWaveVectors(Space _space, int numAtoms, double 
            density, double temperature, String filename, double harmonicFudge,
            int[] compWV, int harmWV){
        super(_space, true);
        
        //Set up some of the joint stuff
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        basis = new BasisMonatomic(space);
        
        integrators = new IntegratorMC[2];
        accumulatorPumps = new DataPump[2];
        meters = new IEtomicaDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];
        
//TARGET    
        // Set up target system - A, 1, hard rod
        PotentialMasterList potentialMasterTarget = new PotentialMasterList(
                this, space);
        boxTarget = new Box(space);
        boxTarget.setNMolecules(species, numAtoms);
        
        Potential2 p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        p2.setBox(boxTarget);
        potentialMasterTarget.addPotential(p2, new IAtomTypeLeaf[] {
                species.getLeafType(), species.getLeafType()});
        
        primitive = new PrimitiveCubic(space, 1.0/density);
        boundaryTarget = new BoundaryRectangularPeriodic(space, numAtoms/density);
        nCells = new int[]{numAtoms};
        boxTarget.setBoundary(boundaryTarget);
        
        CoordinateDefinitionLeaf coordinateDefinitionTarget = new 
                CoordinateDefinitionLeaf(this, boxTarget, primitive, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);
        
        double neighborRange = 1.01/density;
        potentialMasterTarget.setRange(neighborRange);
        // Find neighbors now.  Don't hook up the NieghborListManager since the
        //  neighbors won't change.
        potentialMasterTarget.getNeighborManager(boxTarget).reset();
        
        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget,
                random, temperature);
        integrators[1] = integratorTarget;
        integratorTarget.setBox(boxTarget);
        
        nm = new NormalModes1DHR(space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactoryTarget = nm.getWaveVectorFactory();
        waveVectorFactoryTarget.makeWaveVectors(boxTarget);
        
        changeMove = new MCMoveChangeMode(potentialMasterTarget, random);
        integratorTarget.getMoveManager().addMCMove(changeMove);
        changeMove.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
        changeMove.setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
        
    
    
    
    
    }
    
    
    
    
    public static void main(String args[]){
        

    }
}
