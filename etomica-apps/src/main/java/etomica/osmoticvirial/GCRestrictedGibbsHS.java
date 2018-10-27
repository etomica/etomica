package etomica.osmoticvirial;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2HardSphere;
import etomica.potential.P2Ideal;
import etomica.potential.Potential2;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;

/**
 * Calculate osmotic virial coefficients for Hard-Sphere potential in Grand Canonical ensemble for solvent
 * from Free-Energy Perturbation Approach in restricted grand canonical ensemble
 */
public class GCRestrictedGibbsHS extends Simulation {

    protected IntegratorMC integrator1, integrator2;
    protected IntegratorRGEMC integrator;
    protected MCMoveAtom mcMoveAtom;
    protected MCMoveInsertDelete mcMoveInsertDelete1, mcMoveInsertDelete2;
    protected SpeciesSpheresMono species1;
    protected SpeciesSpheresMono species2;
    protected Box box1, box2;
    protected P2HardSphere potential1, potential12;
    protected Potential2 potential2;
    protected Controller controller;
    protected ActivityIntegrate activityIntegrate;

    /**
     * @param vf reservoir volume fraction of solvent
     * @param q size of solvent divided by size of solute
     */
    public GCRestrictedGibbsHS(double vf, double q, int numAtoms){
        super(Space3D.getInstance());
//        setRandom(new RandomMersenneTwister(1));
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this,space);
        mcMoveInsertDelete1 = new MCMoveInsertDelete(potentialMaster, random, space);
        mcMoveInsertDelete2 = new MCMoveInsertDelete(potentialMaster, random, space);
        species1 = new SpeciesSpheresMono(this, space);
        species2 = new SpeciesSpheresMono(this, space);
        addSpecies(species1);
        addSpecies(species2);

        integrator = new IntegratorRGEMC(random, space, species1);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        double sigma1 = 1; //solute
        double sigma2 = q * sigma1; //solvent
        double sigma12 = (sigma1+sigma2)/2;

        double mu = (8*vf-9*vf*vf+3*vf*vf*vf) / Math.pow((1-vf),3)+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))); //Configurational chemical potential from Carnahan–Starling equation of state
        System.out.println("mu "+ mu+" muig "+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))));

        box1 = new Box(space);
        addBox(box1);
        box1.setBoundary(new BoundaryRectangularPeriodic(space, 3*sigma1));
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
        box2.setBoundary(new BoundaryRectangularPeriodic(space, 3*sigma1));
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
//        potential2 = new P2Ideal(space);
//        System.out.println("AO");
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
        integrator1.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box1).makeMCMoveListener());
        integrator2.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box2).makeMCMoveListener());

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
            params.numSteps = 1000000;
            params.nBlocks = 100;
            params.vf = 0.05;
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

        GCRestrictedGibbsHS sim = new GCRestrictedGibbsHS(vf, q, numAtoms);

        System.out.println("box length "+sim.box1.getBoundary().getBoxSize());
        System.out.println("species1 " +sim.species1.getLeafType().getIndex());
        System.out.println("species2 " +sim.species2.getLeafType().getIndex());

        if (graphics) {
            final String APP_NAME = "SimHard";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box1).canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box1).canvas).setBoundaryFrameColor(Color.BLACK);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box2).canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box2).canvas).setBoundaryFrameColor(Color.BLACK);((DiameterHashByType)simGraphic.getDisplayBox(sim.box1).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box2).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
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

        AccumulatorAverageCovariance acc = new AccumulatorAverageCovariance();
        MCMoveListenerRGE mcMoveListenerRGE = new MCMoveListenerRGE(acc, sim.box1, sim.species1, numAtoms);
        sim.integrator.getMoveEventManager().addListener(mcMoveListenerRGE);
        sim.getController().actionPerformed();

        IData iavg = acc.getData(acc.AVERAGE);
        IData ierr = acc.getData(acc.ERROR);
        IData icor = acc.getData(acc.BLOCK_CORRELATION);
        IData icov = acc.getData(acc.BLOCK_COVARIANCE);

        int n = numAtoms/2+1;

        double[] avg = new double[n];
        double[] err = new double[n];
        double[] cor = new double[n];

        for(int i=0; i<n; i++){
            avg[i] = iavg.getValue(i);
            err[i] = ierr.getValue(i);
            cor[i] = icor.getValue(i);
            System.out.print(String.format("%d avg: %13.6e   err: %11.4e   cor: % 4.2f\n",i, avg[i], err[i], cor[i]));
        }

        for(int i=0; i<n; i++){
            for(int j=i+1; j<n; j++){
                System.out.println("cor: "+icov.getValue(i*n+j)/(Math.sqrt(icov.getValue(i*n+i)*icov.getValue(j*n+j))));
            }
        }

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)*0.001);

     }

    public static class simParams extends ParameterBase{
        public int numAtoms = 3;
        public int numSteps = 200000;
        public int nBlocks = 100;
        public double vf = 0.2;
        public double q = 2.0;
    }
}
