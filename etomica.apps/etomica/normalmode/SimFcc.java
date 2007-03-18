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
import etomica.lattice.BravaisLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.models.hexane.MeterCorrelationMatrix;
import etomica.phase.Phase;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the collective-coordinate S-matrix.
 * No graphic display of simulation. 
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
        phase.setDensity(density);

        if (space.D() == 1) {
            lattice = new LatticeCubicSimple(1,phase.getBoundary().getDimensions().x(0)/numAtoms);
            ((IntegratorHard)integrator).setNullPotential(new P1HardPeriodic(space));
        }
        else {
            lattice = new LatticeCubicFcc();
        }
        config = new ConfigurationLattice(lattice);

        config.initializeCoordinates(phase);

        integrator.setPhase(phase);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        
        //defaults
        int D = 1;
        int nA = 108;
        double density = 1.04;
        if (D == 1) {
            nA = 5;
            density = 0.5;
        }
        String filename = "normal_modes"+D+"D";
        double simTime = 400;
        
        //parse arguments
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
        
        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" hard sphere simulation");
        System.out.println(nA+" atoms at density "+density);
        System.out.println(simTime+" time units");
        System.out.println("output data to "+filename);

        //construct simulation
        SimFcc sim = new SimFcc(Space.getInstance(D), nA, density);
        
        //set up initial configuration and save nominal positions
        Primitive primitive = sim.lattice.getPrimitive();//lattice used to position atoms, not scaled to phase volume
        if (D == 3) {
            primitive = ((LatticeCubicFcc)sim.lattice).getPrimitiveFcc();
        }
        ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice) sim.config
                .getLatticeMemento();//lattice indicating atom positions
        IVector scaling = myLattice.latticeScaling;
        primitive.scaleSize(scaling.x(0));//rescale primitive to indicate actual positioning of atoms

        //set up normal-mode meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(new CoordinateDefinitionLeaf(sim.getSpace()));
        WaveVectorFactory waveVectorFactory;
        if (D == 1) {
            waveVectorFactory = new WaveVectorFactory1D(primitive);
        }
        else if (D == 2) {
            waveVectorFactory = null;
        }
        else {
            waveVectorFactory = new WaveVectorFactoryFcc((PrimitiveFcc)primitive);
        }
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setPhase(sim.phase);

        IntervalActionAdapter fooAdapter = new IntervalActionAdapter(meterNormalMode);
        fooAdapter.setActionInterval(2);
        sim.integrator.addListener(fooAdapter);

        //start simulation
        int nSteps = (int) (simTime / sim.integrator.getTimeStep());
        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();
        
        //normalize averages
        DataGroup normalModeData = (DataGroup)meterNormalMode.getData();
        normalModeData.TE(1.0/(sim.phase.getSpeciesMaster().moleculeCount()*meterNormalMode.getCallCount()));
        int normalDim = meterNormalMode.getCoordinateDefinition().getCoordinateDim();
        
        //write results to file
        IVector[] waveVectors = waveVectorFactory.getWaveVectors();
        double[] coefficients = waveVectorFactory.getCoefficients();
        
        try {
            FileWriter fileWriterQ = new FileWriter(filename+".Q");
            FileWriter fileWriterS = new FileWriter(filename+".S");
            for (int i=0; i<waveVectors.length; i++) {
                fileWriterQ.write(Double.toString(coefficients[i]));
                for (int j=0; j<waveVectors[i].getD(); j++) {
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
    public BravaisLattice lattice;
    public ConfigurationLattice config;
}