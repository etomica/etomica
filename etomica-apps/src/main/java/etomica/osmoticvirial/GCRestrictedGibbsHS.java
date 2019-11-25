package etomica.osmoticvirial;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.cell.PotentialMasterCellMixed;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.HashSet;
import java.util.Set;

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
    protected ActivityIntegrate activityIntegrate;

    /**
     * @param vf reservoir volume fraction of solvent
     * @param q size of solvent divided by size of solute
     */
    public GCRestrictedGibbsHS(double vf, double q, int numAtoms, boolean computeAO, double L, double GMIfac) {
        super(Space3D.getInstance());
//        setRandom(new RandomMersenneTwister(1));
        PotentialMaster potentialMaster;
        potentialMaster = new PotentialMasterCellMixed(this, q, space);
//        potentialMaster = new PotentialMasterCell(this,1,space);
//        potentialMaster = new PotentialMasterMonatomic(this);
        PotentialMasterCell pmc = potentialMaster instanceof PotentialMasterCell ? (PotentialMasterCell) potentialMaster : null;
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

        double mu = Math.log(6*vf/(Math.PI*Math.pow(sigma2,3)));

        potential1 = new P2HardSphere(space, sigma1, false);
        potential12 = new P2HardSphere(space, sigma12, false);
        if (pmc != null) pmc.setCellRange(potentialMaster instanceof PotentialMasterCellMixed ? 2 : 3);
        potentialMaster.setPotentialHard(true);

        if (!computeAO){
            mu+= (8*vf-9*vf*vf+3*vf*vf*vf) / Math.pow((1-vf),3); //Configurational chemical potential from Carnahan–Starling equation of state
            potential2 = new P2HardSphere(space, sigma2, false);
        }
        else{
            potential2 = new P2Ideal(space);
            System.out.println("AO");
        }

        System.out.println("mu "+ mu+" muig "+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))));

        box1 = new Box(space);
        addBox(box1);
        box1.setBoundary(new BoundaryRectangularPeriodic(space, L * sigma1));
        box1.setNMolecules(species1,numAtoms/2);
        integrator1 = new IntegratorMC(this, potentialMaster);
        integrator1.setBox(box1);
        MCMoveManager moveManager = integrator1.getMoveManager();
        if(vf != 0) {
            mcMoveInsertDelete1.setSpecies(species2);
            mcMoveInsertDelete1.setMu(mu);
            moveManager.addMCMove(mcMoveInsertDelete1);
        }
        moveManager.addMCMove(new MCMoveAtom(random, potentialMaster, space));
        integrator.addIntegrator(integrator1);

        box2 = new Box(space);
        addBox(box2);
        box2.setBoundary(new BoundaryRectangularPeriodic(space, L * sigma1));
        box2.setNMolecules(species1,numAtoms - (numAtoms/2));
        integrator2 = new IntegratorMC(this, potentialMaster);
        integrator2.setBox(box2);
        moveManager = integrator2.getMoveManager();
        if(vf != 0) {
            mcMoveInsertDelete2.setSpecies(species2);
            mcMoveInsertDelete2.setMu(mu);
            moveManager.addMCMove(mcMoveInsertDelete2);
        }
        moveManager.addMCMove(new MCMoveAtom(random, potentialMaster, space));
        integrator.addIntegrator(integrator2);
        AtomType leafType1 = species1.getLeafType();
        AtomType leafType2 = species2.getLeafType();

        if (potentialMaster instanceof PotentialMasterCellMixed) {
            ((PotentialMasterCellMixed) potentialMaster).setHandledByCells(species2, true);
            ((PotentialMasterCellMixed) potentialMaster).addUnrangedPotential(leafType1, leafType1, potential1);
            ((PotentialMasterCellMixed) potentialMaster).addUnrangedPotential(leafType1, leafType2, potential12);
            potentialMaster.addPotential(potential2, new AtomType[]{leafType2, leafType2});
        } else {
            potentialMaster.addPotential(potential1, new AtomType[]{leafType1, leafType1});
            potentialMaster.addPotential(potential12, new AtomType[]{leafType1, leafType2});
            potentialMaster.addPotential(potential2, new AtomType[]{leafType2, leafType2});
        }
        integrator.getMCMoveMoleculeExchange().setPotential(leafType1, leafType1, potential1);
        integrator.getMCMoveMoleculeExchange().setPotential(leafType1, leafType2, potential12);
        integrator.getMCMoveMoleculeExchange().setPotential(leafType2, leafType2, potential2);

        // add solvents to the boxes as though there are no solutes (and add solutes as though there
        // are no solvents).  then remove any solvents that overlap a solute
        int nSolvent = (int) (box1.getBoundary().volume() * vf / (q * q * q * Math.PI / 6));
        box1.setNMolecules(species2, nSolvent);
        box2.setNMolecules(species2, nSolvent);

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box1, species1);
        configuration.initializeCoordinates(box1, species2);
        final Set<IMolecule> overlaps = new HashSet<>();
        PotentialCalculation overlapCheck = new PotentialCalculation() {
            @Override
            public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
                double u = potential.energy(atoms);
                if (u < Double.POSITIVE_INFINITY) return;
                IAtom a = atoms.getAtom(0);
                if (a.getType().getSpecies() == species1) a = atoms.getAtom(1);
                overlaps.add(a.getParentGroup());
            }
        };
        IteratorDirective id = new IteratorDirective();
        potentialMaster.calculate(box1, id, overlapCheck);
        for (IMolecule m : overlaps) box1.removeMolecule(m);
        // now rinse and repeat with box2
        configuration.initializeCoordinates(box2, species1);
        configuration.initializeCoordinates(box2, species2);
        overlaps.clear();
        potentialMaster.calculate(box2, id, overlapCheck);
        for (IMolecule m : overlaps) box2.removeMolecule(m);

        // set the global move interval.  GMIfac yields a simulation that spends roughly equal
        // CPU time doing global and not-global moves.  Set GMIfac to achieve desired result.
        double n = (box1.getNMolecules(species2) + box2.getNMolecules(species2)) / 2.0;
        if(n==0){
            integrator.setGlobalMoveInterval(1);
        }
        else {
            integrator.setGlobalMoveInterval(n / GMIfac);
        }

        if (pmc != null) {
            // we have cell listing.  initialize all that.
            integrator.getMoveEventManager().addListener(pmc.getNbrCellManager(box1).makeMCMoveListener());
            integrator.getMoveEventManager().addListener(pmc.getNbrCellManager(box2).makeMCMoveListener());
            integrator1.getMoveEventManager().addListener(pmc.getNbrCellManager(box1).makeMCMoveListener());
            integrator2.getMoveEventManager().addListener(pmc.getNbrCellManager(box2).makeMCMoveListener());
            pmc.getNbrCellManager(box1).assignCellAll();
            pmc.getNbrCellManager(box2).assignCellAll();
        }

    }

    public static void main(String[] args) {
        simParams params = new simParams();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.numAtoms = 2;
            params.numSteps = 5000000000L;
            params.vf = 0.23;
            params.q = 0.2;
            params.computeAO = false;
            params.L = 3;
            params.GMIfac = 3;
        }

        int numAtoms = params.numAtoms;
        long numSteps = params.numSteps;
        double vf = params.vf;
        double q = params.q;
        boolean graphics = false;
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

        GCRestrictedGibbsHS sim = new GCRestrictedGibbsHS(vf, q, numAtoms, computeAO, L, GMIfac);
        int[] seeds = sim.getRandomSeeds();
        System.out.println("Random Seeds:");
        for (int i = 0; i < seeds.length; ++i) {
            System.out.println(seeds[i]);
        }
        System.out.println("species1 " + sim.species1.getLeafType().getIndex());
        System.out.println("species2 " + sim.species2.getLeafType().getIndex());

        if (graphics) {
            final String APP_NAME = "SimHard";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());

            ((DiameterHashByType) simGraphic.getDisplayBox(sim.box1).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            ((DiameterHashByType) simGraphic.getDisplayBox(sim.box2).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            simGraphic.makeAndDisplayFrame(APP_NAME);


            AccumulatorAverageFixed acc = new AccumulatorAverageFixed();
            MCMoveListenerRGE mcMoveListenerRGE = new MCMoveListenerRGE(acc, sim.box1, sim.species1, numAtoms);
            sim.integrator.getMoveEventManager().addListener(mcMoveListenerRGE);
            return;
        }

        sim.activityIntegrate.setMaxSteps(numSteps / 10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);

        double GMI = sim.integrator.getGlobalMoveInterval();
        int targetNumBlocks = 100;
        long blockSize = (long) (numSteps / (targetNumBlocks * GMI));
        if (blockSize == 0) blockSize = 1;
        AccumulatorAverageCovariance acc = new AccumulatorAverageCovariance(blockSize);
        MCMoveListenerRGE mcMoveListenerRGE = new MCMoveListenerRGE(acc, sim.box1, sim.species1, numAtoms);
        sim.integrator.getMoveEventManager().addListener(mcMoveListenerRGE);
        sim.getController().actionPerformed();

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
