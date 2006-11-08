package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.models.hexane.MeterCorrelationMatrix;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * A monatomic fcc hard sphere simulation to test a new energy method.
 * 
 * @author nancycribbin
 * 
 */
public class SimFcc extends Simulation {

    public SimFcc(Space space, int numAtoms, double density) {
        super(space, true, new PotentialMaster(space));

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);

        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);

        integrator = new IntegratorHard(this);

        integrator.setIsothermal(false);
        activityIntegrate = new ActivityIntegrate(this,
                integrator);
        double timeStep = 0.04;
        integrator.setTimeStep(timeStep);
        getController().addAction(activityIntegrate);
        // activityIntegrate.setMaxSteps(nSteps);

        Potential potential = new P2HardSphere(space, defaults.atomSize, false);
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomFactoryMono) species
                .moleculeFactory()).getType();
        potentialMaster.addPotential(potential, new AtomType[] { sphereType,
                sphereType });

        bdry = new BoundaryRectangularPeriodic(this);
        phase.setBoundary(bdry);
        phase.setDensity(density);

        lattice = new LatticeCubicFcc();
        config = new ConfigurationLattice(lattice);
        // config.setRescalingToFitVolume(false);

        config.initializeCoordinates(phase);

        // nan phase.setDensity(1.04);
        integrator.setPhase(phase);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        int nA = 108;
        String filename = "normal_modes110";
        double density = 1.10;
        double simTime = 400;
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            simTime = Double.parseDouble(args[2]);
        }
        if (args.length > 3) {
            nA = Integer.parseInt(args[3]);
        }
        
        System.out.println("Running FCC hard sphere simulation");
        System.out.println(nA+" atoms at density "+density);
        System.out.println(simTime+" time units");
        System.out.println("output data to "+filename);

        SimFcc sim = new SimFcc(Space3D.getInstance(), nA, density);
        
        PrimitiveFcc primitive = sim.lattice.getPrimitiveFcc();
        ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice) sim.config
                .getLatticeMemento();
        Vector scaling = myLattice.latticeScaling;
        primitive.setCubicSize(primitive.getCubicSize()*scaling.x(0));

        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setNormalCoordWrapper(new NormalCoordLeaf(sim.space));
        WaveVectorFactoryFcc waveVectorFactory = new WaveVectorFactoryFcc(primitive);
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setPhase(sim.phase);

        int nSteps = (int) (simTime / sim.integrator.getTimeStep());

        sim.activityIntegrate.setMaxSteps(nSteps);

        IntervalActionAdapter fooAdapter = new IntervalActionAdapter(meterNormalMode);
        fooAdapter.setActionInterval(2);
        sim.integrator.addListener(fooAdapter);
        sim.getController().actionPerformed();
        
        DataGroup normalModeData = (DataGroup)meterNormalMode.getData();
        normalModeData.TE(1.0/(sim.phase.getSpeciesMaster().moleculeCount()*meterNormalMode.getCallCount()));
        int normalDim = meterNormalMode.getNormalCoordWrapper().getNormalDim();
        
        Vector[] waveVectors = meterNormalMode.getWaveVectors();
        
        try {
            FileWriter fileWriterQ = new FileWriter(filename+".Q");
            FileWriter fileWriterS = new FileWriter(filename+".S");
            for (int i=0; i<waveVectors.length; i++) {
                fileWriterQ.write(Double.toString(waveVectors[i].x(0)));
                for (int j=1; j<waveVectors[i].D(); j++) {
                    fileWriterQ.write(" "+waveVectors[i].x(j));
                }
                fileWriterQ.write("\n");
                DataDoubleArray dataS = (DataDoubleArray)normalModeData.getData(i);
                for (int k=0; k<normalDim; k++) {
                    fileWriterS.write(Double.toString(dataS.getValue(k*normalDim)));
                    for (int l=1; l<normalDim; l++) {
                        fileWriterS.write(" "+dataS.getValue(k*normalDim+l));
                    }
                    fileWriterS.write("\n");
                }
            }
            fileWriterQ.close();
            fileWriterS.close();
        }
        catch (IOException e) {
            throw new RuntimeException("Oops, failed to write data "+e);
        }
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMD integrator;
    public ActivityIntegrate activityIntegrate;
    public MeterCorrelationMatrix meterCorrelation;
    public Phase phase;
    public BoundaryRectangularPeriodic bdry;
    public LatticeCubicFcc lattice;
    public ConfigurationLattice config;
}