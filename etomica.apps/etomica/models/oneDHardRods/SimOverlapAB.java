package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
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
import etomica.units.Null;
import etomica.util.ParameterBase;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

public class SimOverlapAB extends Simulation {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimOverlapAB";
    Primitive primitive;
    int[] nCells;
    NormalModes1DHR nm;
    double alpha;       //adjustable parameter
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    ActivityIntegrate activityIntegrate;
    
    
    IntegratorMC[] integrators;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] pumps;
    public DataSource[] meters;
    public IBox boxTarget, boxRef;
    public Boundary boundaryTarget, boundaryRef;
    MCMoveChangeMode changeMove;
    MCMoveConvertMode convertMove;
    MeterPotentialEnergy meterAinB, meterAinA;
    MeterConvertModeBrute meterBinA, meterBinB;
    MeterOverlap meterOverlapInA, meterOverlapInB;
    
    public SimOverlapAB(Space _space, int numAtoms, double density, double 
            temperature, String filename, double harmonicFudge, int awv){
        super(_space, true);
        
        //Set up some of the joint stuff
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        integrators = new IntegratorMC[2];
        pumps = new DataPump[2];
        meters = new DataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];
        
        
        
        //Set up target system   -A
        PotentialMasterList potentialMasterTarget = new 
            PotentialMasterList(this, space);
        boxTarget = new Box(this, space);
        addBox(boxTarget);
        boxTarget.setNMolecules(species, numAtoms);
        
        Potential2 p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        p2.setBox(boxTarget);
        potentialMasterTarget.addPotential(p2, new IAtomTypeLeaf[]
                {species.getLeafType(), species.getLeafType()});
        
        primitive = new PrimitiveCubic(space, 1.0/density);
        boundaryTarget = new BoundaryRectangularPeriodic(space, getRandom(),
                numAtoms/density);
        nCells = new int[]{numAtoms};
        boxTarget.setBoundary(boundaryTarget);
        
        CoordinateDefinitionLeaf coordinateDefinitionTarget = new 
                CoordinateDefinitionLeaf(this, boxTarget, primitive, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);
        
        double neighborRange = 1.01/density;
        potentialMasterTarget.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMasterTarget.getNeighborManager(boxTarget).reset();
        
        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget, 
                random, temperature);
        integrators[1] = integratorTarget;
        integratorTarget.setBox(boxTarget);
        activityIntegrate = new ActivityIntegrate(integratorTarget);
        getController().addAction(activityIntegrate);
        
        nm = new NormalModes1DHR(space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactoryTarget = nm.getWaveVectorFactory();
        waveVectorFactoryTarget.makeWaveVectors(boxTarget);
        
        changeMove = new MCMoveChangeMode(potentialMasterTarget, random);
        integratorTarget.getMoveManager().addMCMove(changeMove);
        changeMove.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
        changeMove.setWaveVectorCoefficients(
                waveVectorFactoryTarget.getCoefficients());
        changeMove.setEigenVectors(nm.getEigenvectors(boxTarget));
        changeMove.setCoordinateDefinition(coordinateDefinitionTarget);
        changeMove.setBox((IBox)boxTarget);
        changeMove.setStepSizeMin(0.001);
        changeMove.setStepSize(0.01);
        
        meterAinA = new MeterPotentialEnergy(potentialMasterTarget);
        meterAinA.setBox(boxTarget);
        
        meterBinA = new MeterConvertModeBrute(potentialMasterTarget, 
                coordinateDefinitionTarget, boxTarget);
        meterBinA.setEigenVectors(nm.getEigenvectors(boxTarget));
        meterBinA.setOmegaSquared(nm.getOmegaSquared(boxTarget));
        meterBinA.setTemperature(temperature);
        meterBinA.setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
        meterBinA.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
        
        meterOverlapInA = new MeterOverlap("MeterOverlapInA", Null.DIMENSION, 
                meterAinA, meterBinA, temperature);
        
        integratorTarget.setBox(boxTarget);
        //nan do we need to reset the potentialmaster neighbormanager here?
        potentialMasterTarget.getNeighborManager(boxTarget).reset();

        
        
        
        
        //Set up reference system - Right now, SystemB - hybrid system
        PotentialMasterList potentialMasterRef = new 
            PotentialMasterList(this, space);
        boxRef = new Box(this, space);
        addBox(boxRef);
        boxRef.setNMolecules(species, numAtoms);
        
        p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        p2.setBox(boxRef);
        potentialMasterRef.addPotential(p2, new IAtomTypeLeaf[]
                {species.getLeafType(), species.getLeafType()});
        
        primitive = new PrimitiveCubic(space, 1.0/density);
        boundaryRef = new BoundaryRectangularPeriodic(space, getRandom(),
                numAtoms/density);
        nCells = new int[]{numAtoms};
        boxRef.setBoundary(boundaryRef);
        
        CoordinateDefinitionLeaf coordinateDefinitionRef = new 
                CoordinateDefinitionLeaf(this, boxRef, primitive, space);
        coordinateDefinitionRef.initializeCoordinates(nCells);
        
        neighborRange = 1.01/density;
        potentialMasterRef.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMasterRef.getNeighborManager(boxRef).reset();
        
        IntegratorMC integratorRef = new IntegratorMC(potentialMasterRef, 
                random, temperature);
        integratorRef.setBox(boxRef);
        activityIntegrate = new ActivityIntegrate(integratorRef);
        getController().addAction(activityIntegrate);
        
        nm = new NormalModes1DHR(space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactoryRef = nm.getWaveVectorFactory();
        waveVectorFactoryRef.makeWaveVectors(boxRef);
        
        
        convertMove = new MCMoveConvertMode(potentialMasterRef, 
                random);
        integratorRef.getMoveManager().addMCMove(convertMove);
        convertMove.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
        convertMove.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
        convertMove.setOmegaSquared(nm.getOmegaSquared(boxRef), 
                waveVectorFactoryRef.getCoefficients());
        convertMove.setEigenVectors(nm.getEigenvectors(boxRef));
        convertMove.setCoordinateDefinition(coordinateDefinitionRef);
        convertMove.setTemperature(temperature);
        convertMove.setBox((IBox)boxRef);
        convertMove.setStepSizeMin(0.001);
        convertMove.setStepSize(0.01);
        
        meterAinB = new MeterPotentialEnergy(potentialMasterRef);
        meterAinB.setBox(boxRef);
       
        meterBinB = new MeterConvertModeBrute(potentialMasterRef,
                coordinateDefinitionRef, boxRef);
        meterBinB.setCoordinateDefinition(coordinateDefinitionRef);
        meterBinB.setEigenVectors(nm.getEigenvectors(boxRef));
        meterBinB.setOmegaSquared(nm.getOmegaSquared(boxRef));
        meterBinB.setTemperature(temperature);
        meterBinB.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
        meterBinB.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
        
        
        meterOverlapInB = new MeterOverlap("MeterOverlapInB", Null.DIMENSION, 
                meterAinB, meterBinB, temperature);
        
        integratorRef.setBox(boxRef);
        potentialMasterRef.getNeighborManager(boxRef).reset();
        
        
        
        //Set up the rest of the joint stuff
        setAffectedWaveVector(awv);
        
        integratorOverlap = new IntegratorOverlap(random, new 
                IntegratorMC[]{integratorRef, integratorTarget});
        
        
        
        
    }
    
    
    
    public static void main(String args[]){
        
        
    }
    
    
    public void setAffectedWaveVector(int awv){
        convertMove.setConvertedWaveVector(awv);
        meterBinA.setConvertedWV(awv);
        meterBinB.setConvertedWV(awv);
    }
    public static class SimOverlapABParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 0.5;
        public int D = 1;
        public long numSteps = 1000;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public double temperature = 1.0;
        public int affectedWV = 16;
    }
    
}
