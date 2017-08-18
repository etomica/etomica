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
import etomica.data.meter.MeterWidomInsertion;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Implements Free-Energy Perturbation Approach (Widom's Insertion method) for calculation of osmotic virial coefficient
 * for Hard-Sphere potential. GC ensemble
 * Created by aksharag on 6/16/17.
 */
public class fepHSgc extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public MCMoveInsertDelete mcMoveInsertDelete ;
    public SpeciesSpheresMono species1;
    public SpeciesSpheresMono species2;
    public Box box;
    public P2HardSphere potential1, potential2, potential12;
    public Controller controller;
    public ActivityIntegrate activityIntegrate;

    public fepHSgc(int numAtoms, double vf, double q, boolean computez2z1, boolean computez3z2){
        super(Space3D.getInstance());
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this,space);

        integrator = new IntegratorMC(this, potentialMaster);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        mcMoveInsertDelete = new MCMoveInsertDelete(potentialMaster, random, space);

        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getMoveManager().addMCMove(mcMoveInsertDelete);

        double sigma1 = 1.0; //solute
        double sigma2 = q * sigma1; //solvent
        double sigma12 = (sigma1+sigma2)/2;

        species1 = new SpeciesSpheresMono(this, space);
        species2 = new SpeciesSpheresMono(this, space);
        addSpecies(species1);
        addSpecies(species2);
        box = new Box(space);
        addBox(box);
        box.setBoundary(new BoundaryRectangularPeriodic(space, 4*sigma1));

        mcMoveInsertDelete.setSpecies(species2);
        double mu = (8*vf-9*vf*vf+3*vf*vf*vf) / Math.pow((1-vf),3)+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))); //Configurational chemical potential from Carnahan–Starling equation of state
        System.out.println("mu "+ mu+" muig "+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))));
        mcMoveInsertDelete.setMu(mu);

        if (computez2z1){box.setNMolecules(species1,1);}
        else if (computez3z2){box.setNMolecules(species1,2);}

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

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);

        integrator.setBox(box);
        potentialMaster.getNbrCellManager(box).assignCellAll();

    }

    public static void main(String[] args){
        simParams params = new simParams();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numAtoms = 3;
            params.numSteps = 5000;
            params.nBlocks = 1000;
            params.vf = 0.1;
            params.q = 0.2;
            params.computez2z1 = true;
            params.computez3z2 = false;
        }

        int numAtoms = params.numAtoms;
        int numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double vf = params.vf;
        double q = params.q;
        boolean computez2z1 = params.computez2z1;
        boolean computez3z2 = params.computez3z2;
        boolean graphics = false;

        long numSamples = numSteps/numAtoms;
        long samplesPerBlock = numSamples/nBlocks;
        if (samplesPerBlock == 0) samplesPerBlock = 1;

        System.out.println("Hard Sphere OV FEP GC");

        if(computez2z1){
            System.out.println("**z2_z1**");
        }
        else if(computez3z2){
            System.out.println("**z3_z2**");
        }
        else{
            System.out.println("**z1_z0**");
        }

        System.out.println(numAtoms+" atoms, "+numSteps+" steps");
        System.out.println("vol fraction: "+vf);
        System.out.println("q: "+q);
        System.out.println(nBlocks+" blocks");

        long t1 = System.currentTimeMillis();

        fepHSgc sim = new fepHSgc(numAtoms, vf, q, computez2z1, computez3z2);

        System.out.println("box length "+sim.box.getBoundary().getBoxSize());

        if (graphics) {
            final String APP_NAME = "SimHard";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());

            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            simGraphic.makeAndDisplayFrame(APP_NAME);

            MeterWidomInsertion meterinsert = new MeterWidomInsertion(sim.space,sim.getRandom());
            //meterinsert.setNInsert(50);
            meterinsert.setSpecies(sim.species1);
            meterinsert.setIntegrator(sim.integrator);

            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(samplesPerBlock);
            DataPumpListener pump = new DataPumpListener(meterinsert, acc, numAtoms);
            sim.integrator.getEventManager().addListener(pump);

            return;
        }

        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);

        MeterWidomInsertion meterinsert = new MeterWidomInsertion(sim.space,sim.getRandom());
        //meterinsert.setNInsert(50);
        meterinsert.setSpecies(sim.species1);
        meterinsert.setIntegrator(sim.integrator);

        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pump = new DataPumpListener(meterinsert, acc, numAtoms);
        sim.integrator.getEventManager().addListener(pump);

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
        public int numAtoms = 500;
        public int numSteps = 200000;
        public int nBlocks = 1000;
        public double vf = 0.2;
        public double q = 2.0;
        public boolean computez2z1 = false;
        public boolean computez3z2 = false;

    }
}
