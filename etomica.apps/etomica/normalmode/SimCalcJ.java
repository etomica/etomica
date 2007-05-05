package etomica.normalmode;

import Jama.Matrix;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationLatticeSimple;
import etomica.integrator.IntegratorMD;
import etomica.lattice.BravaisLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation class of hard spheres in 1D or 3D that calculates the dq/dx
 * Jacobian.
 */
public class SimCalcJ extends Simulation {

    public SimCalcJ(Space space, int numAtoms) {
        super(space, true, new PotentialMaster(space));

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);

        bdry = new BoundaryRectangularPeriodic(space, getRandom(), 1);
        phase = new Phase(bdry);
        addPhase(phase);
        phase.getAgent(species).setNMolecules(numAtoms);

        IVector dimensions = bdry.getDimensions();
        if (space.D() == 1) {
            dimensions.E(numAtoms);
        }
        else {
            dimensions.E(Math.round(Math.pow(numAtoms/4, 1.0/3.0)));
        }
        bdry.setDimensions(dimensions);
        phase.setBoundary(bdry);

        if (space.D() == 1) {
            lattice = new LatticeCubicSimple(1, phase.getBoundary().getDimensions().x(0) / numAtoms);
        } else {
            lattice = new LatticeCubicFcc(1);
        }
        config = new ConfigurationLatticeSimple(lattice);

        config.initializeCoordinates(phase);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 1;
        int nA = 4;

        // parse arguments
        if (args.length > 0) {
            nA = Integer.parseInt(args[0]);
        }

        System.out.println("Calculating "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere Jacobian for N="+nA);

        // construct simulation
        SimCalcJ sim = new SimCalcJ(Space.getInstance(D), nA);

        // set up initial configuration and save nominal positions
        Primitive primitive = sim.lattice.getPrimitive();// lattice used to position atoms, not scaled to phase volume
        if (D == 3) {
            primitive = ((LatticeCubicFcc) sim.lattice).getPrimitiveFcc();
        }

        // set up normal-mode meter
        WaveVectorFactory waveVectorFactory;
        if (D == 1) {
            waveVectorFactory = new WaveVectorFactory1D();
        } else if (D == 2) {
            waveVectorFactory = null;
        } else {
            waveVectorFactory = new WaveVectorFactoryFcc((PrimitiveFcc) primitive);
        }
        CalcJacobian meterJacobian = new CalcJacobian(D);
        meterJacobian.setWaveVectorFactory(waveVectorFactory);
        meterJacobian.setPhase(sim.phase);

        double[][] jac = meterJacobian.getJacobian();
        
        for (int i=0; i<jac.length; i++) {
            for (int j=0; j<jac[0].length; j++) {
                System.out.print(jac[i][j]+"  ");
            }
            System.out.println();
        }
        
        Matrix m = new Matrix(jac);
        double d = m.det();
        double d2 = Math.abs(d) / (Math.pow(nA,(D*(nA-1))/2.0));
        System.out.println("raw det = "+d+" (log2 det = "+Math.log(Math.abs(d))/Math.log(2)+")");
        System.out.println("dx/dq det = "+1/d2+" (log2 det = "+Math.log(1.0/d2)/Math.log(2)+")");
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMD integrator;
    public ActivityIntegrate activityIntegrate;
    public Phase phase;
    public BoundaryRectangularPeriodic bdry;
    public BravaisLattice lattice;
    public ConfigurationLatticeSimple config;
}