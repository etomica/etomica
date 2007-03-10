package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.config.ConfigurationLattice;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
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
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;

/**
 * A monatomic fcc hard sphere simulation to test a new energy method.
 * 
 * @author nancycribbin
 * 
 */
public class TestFcc extends Simulation {

    public TestFcc(Space space, int numAtoms) {
        super(space, true, new PotentialMaster(space));

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;
        defaults.ignoreOverlap = true;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);
        getSpeciesRoot().addSpecies(species);

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
        phase.setDensity(1.04);

        lattice = new LatticeCubicFcc();
        config = new ConfigurationLattice(lattice);
        // config.setRescalingToFitVolume(false);

        config.initializeCoordinates(phase);

        // nan phase.setDensity(1.04);
        integrator.setPhase(phase);


        System.out.println(phase.getDensity());
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        int nA = 108;
        boolean graphic = false;
        TestFcc sim = new TestFcc(Space3D.getInstance(), nA);
        
        if(graphic){
            SimulationGraphic simG = new SimulationGraphic(sim);
            simG.makeAndDisplayFrame();
        } else {
            // nan this section is a patch
            // first we find out the scaling used in
            // ConfigurationLattice/LatticeCubicFcc
            // then, we create a primitive fcc lattice, and scale it so we can use
            // it in pri.
            PrimitiveFcc primitive = sim.lattice.getPrimitiveFcc();
            ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice) sim.config
                    .getLatticeMemento();
            IVector scaling = myLattice.latticeScaling;
            primitive.setCubicSize(primitive.getCubicSize()*scaling.x(0));

            MeterNormalMode meterNormalMode = new MeterNormalMode();
            meterNormalMode.setWaveVectorFactory(new WaveVectorFactoryFcc(primitive));
            meterNormalMode.setNormalCoordWrapper(new NormalCoordLeaf(sim.getSpace()));
            meterNormalMode.setPhase(sim.phase);

            if (false) {
                // set up a contrived wave
                sim.phase.setDimensions(new Vector3D(6,6,6));
                ConfigurationLattice config = new ConfigurationLattice(sim.lattice);
                // config.setRescalingToFitVolume(false);
    
                config.initializeCoordinates(sim.phase);
    
                // nan this section is a patch
                // first we find out the scaling used in
                // ConfigurationLattice/LatticeCubicFcc
                // then, we create a primitive fcc lattice, and scale it so we can use
                // it in pri.
                primitive = sim.lattice.getPrimitiveFcc();
                myLattice = (ConfigurationLattice.MyLattice) config
                        .getLatticeMemento();
                scaling = myLattice.latticeScaling;
                primitive.setCubicSize(primitive.getCubicSize()*scaling.x(0));
                System.out.println(primitive.getCubicSize());
                if (primitive.getCubicSize() > Math.sqrt(2)+0.0001 || primitive.getCubicSize() < Math.sqrt(2)-0.0001) {
                    throw new RuntimeException("It should have been 2");
                }
                meterNormalMode.setWaveVectorFactory(new WaveVectorFactoryFcc(primitive));
                AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
                iterator.setPhase(sim.phase);
                iterator.reset();
                while (iterator.hasNext()) {
                    AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
                    IVector pos = atom.getCoord().getPosition();
                    pos.setX(0, pos.x(0)-0.5);
                }

                meterNormalMode.setPhase(sim.phase);
                
                iterator.setPhase(sim.phase);
                iterator.reset();
                while (iterator.hasNext()) {
                    AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
                    IVector pos = atom.getCoord().getPosition();
                    if (Math.round(pos.x(0)+0.5) % 2 == 0) {
                        pos.setX(1,pos.x(1)+0.001);
                    }
                    else {
                        pos.setX(1,pos.x(1)-0.001);
                    }
                }
                meterNormalMode.actionPerformed();
            }            
            
            double simTime = 400.0;
            int nSteps = (int) (simTime / sim.integrator.getTimeStep());

            String filename = "normal_modes4000";
            if (args.length > 0) {
                filename = args[0];
            }
            
            sim.activityIntegrate.setMaxSteps(nSteps);

            IntervalActionAdapter fooAdapter = new IntervalActionAdapter(meterNormalMode);
            fooAdapter.setActionInterval(2);
            sim.integrator.addListener(fooAdapter);
            sim.getController().actionPerformed();
            
            DataGroup normalModeData = (DataGroup)meterNormalMode.getData();
            normalModeData.TE(1.0/(sim.phase.getSpeciesMaster().moleculeCount()*meterNormalMode.getCallCount()));
            int normalDim = meterNormalMode.getNormalCoordWrapper().getNormalDim();
            
            IVector[] waveVectors = meterNormalMode.getWaveVectors();
            
            try {
                FileWriter fileWriterQ = new FileWriter(filename+".Q");
                FileWriter fileWriterS = new FileWriter(filename+".S");
                for (int i=0; i<waveVectors.length; i++) {
                    fileWriterQ.write(Double.toString(waveVectors[i].x(0)));
                    for (int j=1; j<waveVectors[i].getD(); j++) {
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
        
        System.out.println("Peace be unto you.");

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