package etomica.osmoticvirial;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.integrator.mcmove.MCMoveInsertDeleteFasterer;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.potential.*;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.NeighborManagerSimple;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.HashSet;
import java.util.Set;

/**
 * Calculate osmotic virial coefficients for Hard-Sphere potential in Grand Canonical ensemble for solvent
 * from Free-Energy Perturbation Approach in restricted grand canonical ensemble
 */
public class GCRestrictedGibbsHSFasterer extends Simulation {

    protected IntegratorMCFasterer integrator1, integrator2;
    protected IntegratorRGEMCFasterer integrator;
    protected MCMoveAtomFasterer mcMoveAtom;
    protected MCMoveInsertDeleteFasterer mcMoveInsertDelete1, mcMoveInsertDelete2;
    protected SpeciesGeneral species1;
    protected SpeciesGeneral species2;
    protected Box box1, box2;
    protected P2HardGeneric potential1, potential12;
    protected Potential2Soft potential2;


    /**
     * @param vf reservoir volume fraction of solvent
     * @param q size of solvent divided by size of solute
     */
    public GCRestrictedGibbsHSFasterer(double vf, double q, int numAtoms, boolean computeAO, double L, double GMIfac) {
        super(Space3D.getInstance());
//        setRandom(new RandomMersenneTwister(1));

        species1 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species1);
        addSpecies(species2);

        double sigma1 = 1; //solute
        double sigma2 = q * sigma1; //solvent
        double sigma12 = (sigma1+sigma2)/2;

        box1 = new Box(new BoundaryRectangularPeriodic(space, L * sigma1), space);
        addBox(box1);
        box1.setNMolecules(species1, numAtoms / 2);
        box2 = new Box(new BoundaryRectangularPeriodic(space, L * sigma1), space);
        addBox(box2);
        box2.setNMolecules(species1, numAtoms - (numAtoms / 2));

        potential1 = P2HardSphere.makePotential(sigma1);
        potential12 = P2HardSphere.makePotential(sigma12);

        double mu = Math.log(6*vf/(Math.PI*Math.pow(sigma2,3)));

        if (!computeAO){
            mu+= (8*vf-9*vf*vf+3*vf*vf*vf) / Math.pow((1-vf),3); //Configurational chemical potential from Carnahan–Starling equation of state
            potential2 = P2HardSphere.makePotential(sigma2);
        }
        else{
            potential2 = new P2Ideal(space);
            System.out.println("AO");
        }

        System.out.println("mu "+ mu+" muig "+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))));

        PotentialComputePair potentialMaster1, potentialMaster2;
        NeighborManager neighborManager1, neighborManager2;
        if (vf != 0) {
            AtomType solvent = species2.getLeafType();
            AtomType solute = species1.getLeafType();
            neighborManager1 = new NeighborManagerCellMixed(getSpeciesManager(), box1, 2, BondingInfo.noBonding(), solvent);
            neighborManager2 = new NeighborManagerCellMixed(getSpeciesManager(), box2, 2, BondingInfo.noBonding(), solvent);
            potentialMaster1 = new PotentialComputePair(getSpeciesManager(), box1, neighborManager1);
            potentialMaster1.setPairPotential(solvent, solvent, potential2);
            potentialMaster1.setPairPotential(solute, solvent, potential12);
        }
        else{
            neighborManager1 = new NeighborManagerSimple(box1);
            neighborManager2 = new NeighborManagerSimple(box2);
            potentialMaster1 = new PotentialComputePair(getSpeciesManager(), box1, neighborManager1);
        }
        potentialMaster1.setPairPotential(species1.getLeafType(), species1.getLeafType(), potential1);
        potentialMaster2 = new PotentialComputePair(getSpeciesManager(), box2, neighborManager2, potentialMaster1.getPairPotentials());

        mcMoveInsertDelete1 = new MCMoveInsertDeleteFasterer(potentialMaster1, random, space);
        mcMoveInsertDelete2 = new MCMoveInsertDeleteFasterer(potentialMaster2, random, space);

        integrator = new IntegratorRGEMCFasterer(random, space, species1, potentialMaster1.getPairPotentials());
        this.getController().addActivity(new ActivityIntegrate(integrator));

        integrator1 = new IntegratorMCFasterer(potentialMaster1, random, 1, box1);
        MCMoveManager moveManager = integrator1.getMoveManager();
        if(vf != 0) {
            mcMoveInsertDelete1.setSpecies(species2);
            mcMoveInsertDelete1.setMu(mu);
//            moveManager.addMCMove(mcMoveInsertDelete1);
        }
        moveManager.addMCMove(new MCMoveAtomFasterer(random, potentialMaster1, box1));
        integrator.addIntegrator(integrator1, neighborManager1);

        integrator2 = new IntegratorMCFasterer(potentialMaster2, random, 1, box2);
        moveManager = integrator2.getMoveManager();
        if(vf != 0) {
            mcMoveInsertDelete2.setSpecies(species2);
            mcMoveInsertDelete2.setMu(mu);
//            moveManager.addMCMove(mcMoveInsertDelete2);
        }
        moveManager.addMCMove(new MCMoveAtomFasterer(random, potentialMaster2, box2));
        integrator.addIntegrator(integrator2, neighborManager2);
        AtomType leafType1 = species1.getLeafType();
        AtomType leafType2 = species2.getLeafType();

        // add solvents to the boxes as though there are no solutes (and add solutes as though there
        // are no solvents).  then remove any solvents that overlap a solute
        int nSolvent = (int) (box1.getBoundary().volume() * vf / (q * q * q * Math.PI / 6));
        box1.setNMolecules(species2, nSolvent);
        box2.setNMolecules(species2, nSolvent);

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box1, species1);
        configuration.initializeCoordinates(box1, species2);

        potentialMaster1.init();
        potentialMaster2.init();

        if (vf != 0) {
            final Set<IMolecule> overlaps = new HashSet<>();
            PotentialCallback overlapCheck = new PotentialCallback() {
                @Override
                public void pairCompute(int i, int j, Vector dr, double[] u012) {
                    if (u012[0] < Double.POSITIVE_INFINITY) return;
                    IAtom a = box1.getLeafList().get(i);
                    if (a.getType().getSpecies() == species1) a = box1.getLeafList().get(j);
                    overlaps.add(a.getParentGroup());
                }
            };

            potentialMaster1.computeAll(false, overlapCheck);
            for (IMolecule m : overlaps) box1.removeMolecule(m);

            overlaps.clear();
            overlapCheck = new PotentialCallback() {
                @Override
                public void pairCompute(int i, int j, Vector dr, double[] u012) {
                    if (u012[0] < Double.POSITIVE_INFINITY) return;
                    IAtom a = box2.getLeafList().get(i);
                    if (a.getType().getSpecies() == species1) a = box2.getLeafList().get(j);
                    overlaps.add(a.getParentGroup());
                }
            };

            potentialMaster2.computeAll(false, overlapCheck);
            for (IMolecule m : overlaps) box2.removeMolecule(m);
        }

        // set the global move interval.  GMIfac yields a simulation that spends roughly equal
        // CPU time doing global and not-global moves.  Set GMIfac to achieve desired result.
        double n = (box1.getNMolecules(species2) + box2.getNMolecules(species2)) / 2.0;
        if(n==0){
            integrator.setGlobalMoveInterval(1);
        }
        else {
            integrator.setGlobalMoveInterval(n / GMIfac);
        }
    }

    public static void main(String[] args) {
        simParams params = new simParams();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.numAtoms = 2;
            params.numSteps = 5000000000L;
            params.vf = 0.05;
            params.q = 0.2;
            params.computeAO = false;
            params.L = 3;
            params.GMIfac = 3;
        }

        int numAtoms = params.numAtoms;
        long numSteps = params.numSteps;
        double vf = params.vf;
        double q = params.q;
        boolean graphics = true;
        boolean computeAO = params.computeAO;
        double L = params.L;
        double GMIfac = params.GMIfac;

        System.out.println("Hard Sphere OV FEP GC");
        System.out.println(numSteps + " steps long");
        System.out.println("vol fraction: " + vf);
        System.out.println("q: " + q);
        System.out.println("total no of solutes " + numAtoms);
        System.out.println("box length " + L);
        System.out.println("GMI factor " + GMIfac);

        long t1 = System.currentTimeMillis();

        GCRestrictedGibbsHSFasterer sim = new GCRestrictedGibbsHSFasterer(vf, q, numAtoms, computeAO, L, GMIfac);
        int[] seeds = sim.getRandomSeeds();
        System.out.print("Random Seeds:  ");
        for (int i = 0; i < seeds.length; ++i) {
            System.out.print(" "+seeds[i]);
        }
        System.out.println();
        System.out.println("species1 " + sim.species1.getLeafType().getIndex());
        System.out.println("species2 " + sim.species2.getLeafType().getIndex());

        if (graphics) {
            final String APP_NAME = "SimHard";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);

            ((DiameterHashByType) simGraphic.getDisplayBox(sim.box1).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            ((DiameterHashByType) simGraphic.getDisplayBox(sim.box2).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            simGraphic.makeAndDisplayFrame(APP_NAME);


            AccumulatorAverageFixed acc = new AccumulatorAverageFixed();
            MCMoveListenerRGE mcMoveListenerRGE = new MCMoveListenerRGE(acc, sim.box1, sim.species1, numAtoms);
            sim.integrator.getMoveEventManager().addListener(mcMoveListenerRGE);
            return;
        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps / 10));

        sim.integrator.getMoveManager().setEquilibrating(false);

        double GMI = sim.integrator.getGlobalMoveInterval();
        int targetNumBlocks = 100;
        long blockSize = (long) (numSteps / (targetNumBlocks * GMI));
        if (blockSize == 0) blockSize = 1;
        AccumulatorAverageCovariance acc = new AccumulatorAverageCovariance(blockSize);
        MCMoveListenerRGE mcMoveListenerRGE = new MCMoveListenerRGE(acc, sim.box1, sim.species1, numAtoms);
        sim.integrator.getMoveEventManager().addListener(mcMoveListenerRGE);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        System.out.println("block count " + acc.getBlockCount());
        IData iavg = acc.getData(acc.AVERAGE);
        IData ierr = acc.getData(acc.ERROR);
        IData icor = acc.getData(acc.BLOCK_CORRELATION);
        IData icov = acc.getData(acc.BLOCK_COVARIANCE);

        int n = numAtoms / 2 + 1;

        double[] avg = new double[n];
        double[] err = new double[n];
        double[] cor = new double[n];

        for (int i = 0; i < n; i++) {
            avg[i] = iavg.getValue(i);
            err[i] = ierr.getValue(i);
            cor[i] = icor.getValue(i);
            System.out.print(String.format("%d avg: %13.6e   err: %11.4e   cor: % 4.2f\n", i, avg[i], err[i], cor[i]));
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.println("corr " + i + " " + j + " " + icov.getValue(i * n + j));
            }
        }
        if (false){
            int numSwaps[][] = sim.integrator.getMCMoveMoleculeExchange().numSwaps;
            for (int i = 0; i < numSwaps.length; i++) {
                if (numSwaps[i][0] > 0)
                    System.out.println(i + " " + numSwaps[i][0] + " " + numSwaps[i][1]);
            }
        }

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)*0.001);
     }

    public static class simParams extends ParameterBase {
        public int numAtoms = 3;
        public long numSteps = 5000000000L;
        public double vf = 0.2;
        public double q = 0.2;
        public boolean computeAO = false;
        public double L = 3.0;
        public double GMIfac = 3;
    }
}
