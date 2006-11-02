package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import Jama.Matrix;

import etomica.action.PhaseImposePbc;
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
import etomica.integrator.IntegratorMC;
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
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation to sample harmonic potential
 */
public class TestHarmonic extends Simulation {

    public TestHarmonic(Space space, int numAtoms, String filename) {
        super(space, true, new PotentialMaster(space));

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);

        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);

        integrator = new IntegratorMC(this);

        activityIntegrate = new ActivityIntegrate(this,
                integrator);
        getController().addAction(activityIntegrate);

        MCMoveHarmonic move = new MCMoveHarmonic(potentialMaster);
        integrator.getMoveManager().addMCMove(move);
        
        bdry = new BoundaryRectangularPeriodic(this);
        phase.setBoundary(bdry);
        phase.setDensity(1.2);

        PhaseImposePbc makeperiodic = new PhaseImposePbc(phase);
        integrator.addListener(makeperiodic);

        double[][] eigenValues = ArrayReader1D.getFromFile(filename+".val");
        Vector[] q = ArrayReader1D.getVectorsFromFile(filename+".Q");
        double[][][] eigenvectors = ArrayReader2D.getFromFile(filename+".vec");
        
        move.setEigenValues(eigenValues);
        move.setEigenVectors(eigenvectors);
        move.setWaveVectors(q);
        
        lattice = new LatticeCubicFcc();
        ConfigurationLattice config = new ConfigurationLattice(lattice);
        // config.setRescalingToFitVolume(false);

        config.initializeCoordinates(phase);

        move.setPhase(phase);
        
        integrator.setPhase(phase);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        int nA = 108;
        boolean graphic = true;
        String filename = "normal_modes_d12";
        if (args.length > 0) {
            filename = args[0];
        }
        TestHarmonic sim = new TestHarmonic(Space3D.getInstance(), nA, filename);
        
        if(graphic){
            SimulationGraphic simG = new SimulationGraphic(sim);
            simG.makeAndDisplayFrame();
        } else {

            int nSteps = 100000;

            sim.activityIntegrate.setMaxSteps(nSteps);

        }
        
        System.out.println("Peace be unto you.");

    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Phase phase;
    public BoundaryRectangularPeriodic bdry;
    public LatticeCubicFcc lattice;
    public PrimitiveFcc primitive;
}