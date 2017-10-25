package etomica.osmoticvirial;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

/**
 * Calculate osmotic virial coefficients for Hard-Sphere potential in Grand Canonical ensemble for solvent
 * from Free-Energy Perturbation Approach in restricted grand canonical ensemble
 */
public class fepHSgcRGE extends Simulation {

    public IntegratorMC integrator1, integrator2;
    public IntegratorRGEMC integrator;
    public MCMoveAtom mcMoveAtom;
    public MCMoveInsertDelete mcMoveInsertDelete1, mcMoveInsertDelete2;
    public SpeciesSpheresMono species1;
    public SpeciesSpheresMono species2;
    public Box box1, box2;
    public P2HardSphere potential1, potential2, potential12;
    public Controller controller;
    public ActivityIntegrate activityIntegrate;

    /**
     * @param vf reservoir volume fraction of solvent
     * @param q size of solvent divided by size of solute
     */
    public fepHSgcRGE(double vf, double q, int numAtoms){
        super(Space3D.getInstance());
        setRandom(new RandomMersenneTwister(1));
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this,space);
        mcMoveInsertDelete1 = new MCMoveInsertDelete(potentialMaster, random, space);
        mcMoveInsertDelete2 = new MCMoveInsertDelete(potentialMaster, random, space);
        integrator = new IntegratorRGEMC(random,space);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        double sigma1 = 1.0; //solute
        double sigma2 = q * sigma1; //solvent
        double sigma12 = (sigma1+sigma2)/2;

        double mu = (8*vf-9*vf*vf+3*vf*vf*vf) / Math.pow((1-vf),3)+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))); //Configurational chemical potential from Carnahan–Starling equation of state
        System.out.println("mu "+ mu+" muig "+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))));

        species1 = new SpeciesSpheresMono(this, space);
        species2 = new SpeciesSpheresMono(this, space);
        addSpecies(species1);
        addSpecies(species2);

        box1 = new Box(space);
        addBox(box1);
        box1.setBoundary(new BoundaryRectangularPeriodic(space, 4*sigma1));
        box1.setNMolecules(species1,numAtoms/2);
        integrator1 = new IntegratorMC(this, potentialMaster);
        integrator1.setBox(box1);
        MCMoveManager moveManager = integrator1.getMoveManager();
        mcMoveInsertDelete1.setSpecies(species2);
        mcMoveInsertDelete1.setMu(mu);
        moveManager.addMCMove(mcMoveInsertDelete1);
        moveManager.addMCMove(new MCMoveAtom(random, potentialMaster, space));
        integrator.addIntegrator(integrator1);

        box2 = new Box(space);
        addBox(box2);
        box2.setBoundary(new BoundaryRectangularPeriodic(space, 4*sigma1));
        box2.setNMolecules(species1,numAtoms - (numAtoms/2));
        integrator2 = new IntegratorMC(this, potentialMaster);
        integrator2.setBox(box2);
        moveManager = integrator2.getMoveManager();
        mcMoveInsertDelete2.setSpecies(species2);
        mcMoveInsertDelete2.setMu(mu);
        moveManager.addMCMove(mcMoveInsertDelete2);
        moveManager.addMCMove(new MCMoveAtom(random, potentialMaster, space));
        integrator.addIntegrator(integrator2);

        potential1 = new P2HardSphere(space, sigma1, false);
        potential2 = new P2HardSphere(space, sigma2, false);
        potential12 = new P2HardSphere(space, sigma12, false);
        potentialMaster.setCellRange(3);
        potentialMaster.setRange(potential1.getRange());

        AtomType leafType1 = species1.getLeafType();
        AtomType leafType2 = species2.getLeafType();

        potentialMaster.addPotential(potential1, new AtomType[]{leafType1, leafType1});
        potentialMaster.addPotential(potential12, new AtomType[]{leafType1, leafType2});
        potentialMaster.addPotential(potential2, new AtomType[]{leafType2, leafType2});

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box1).makeMCMoveListener());
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box2).makeMCMoveListener());

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box1);
        configuration.initializeCoordinates(box2);

        potentialMaster.getNbrCellManager(box1).assignCellAll();
        potentialMaster.getNbrCellManager(box2).assignCellAll();

    }

    public static void main(String[] args){
        simParams params = new simParams();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numAtoms = 2;
            params.numSteps = 50;
            params.nBlocks = 1000;
            params.vf = 0.1;
            params.q = 0.2;
        }

        int numAtoms = params.numAtoms;
        int numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double vf = params.vf;
        double q = params.q;
        boolean graphics = true;

        long numSamples = numSteps/3;
        long samplesPerBlock = numSamples/nBlocks;
        if (samplesPerBlock == 0) samplesPerBlock = 1;

        System.out.println("Hard Sphere OV FEP GC");
        System.out.println(numSteps+" steps");
        System.out.println("vol fraction: "+vf);
        System.out.println("q: "+q);
        System.out.println(nBlocks+" blocks");
        System.out.println("total no of solutes"+ numAtoms);

        long t1 = System.currentTimeMillis();

        fepHSgcRGE sim = new fepHSgcRGE(vf, q, numAtoms);

        System.out.println("box length "+sim.box1.getBoundary().getBoxSize());

        if (graphics) {
            final String APP_NAME = "SimHard";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());

            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box1).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            simGraphic.makeAndDisplayFrame(APP_NAME);


            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(samplesPerBlock);
            MCMoveListenerRGE mcMoveListenerRGE = new MCMoveListenerRGE(acc, sim.box1, sim.species1, numAtoms);
            sim.integrator.getMoveEventManager().addListener(mcMoveListenerRGE);
            return;
        }

        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);

        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(samplesPerBlock);
        MCMoveListenerRGE mcMoveListenerRGE = new MCMoveListenerRGE(acc, sim.box1, sim.species1, numAtoms);
        sim.integrator.getMoveEventManager().addListener(mcMoveListenerRGE);
        sim.getController().actionPerformed();

        IData iavg = acc.getData(AccumulatorAverage.AVERAGE);
        IData ierr = acc.getData(AccumulatorAverage.ERROR);
        IData icor = acc.getData(AccumulatorAverage.BLOCK_CORRELATION);

        double avg = iavg.getValue(0);
        double err = ierr.getValue(0);
        double cor = icor.getValue(0);

        System.out.print(String.format("avg: %13.6e   err: %11.4e   cor: % 4.2f\n", avg, err, cor));
        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)*0.001);

     }

    public static class simParams extends ParameterBase{
        public int numAtoms = 3;
        public int numSteps = 200000;
        public int nBlocks = 1000;
        public double vf = 0.2;
        public double q = 2.0;

    }
}
