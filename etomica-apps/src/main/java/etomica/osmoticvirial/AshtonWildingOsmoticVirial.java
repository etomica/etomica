package etomica.osmoticvirial;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPumpListener;
import etomica.data.histogram.HistogramSimple;
import etomica.data.meter.MeterNMolecules;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.DoubleRange;
import etomica.molecule.IMolecule;
import etomica.potential.BondingInfo;
import etomica.potential.IPotential2;
import etomica.potential.P2HardSphere;
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
import etomica.util.random.RandomMersenneTwister;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * Implements Ashton and Wilding method for calculation of osmotic virial coefficient (B2, B3)
 * for Hard-Sphere model as described in the paper DOI: 10.1063/1.4883718.
 */
public class AshtonWildingOsmoticVirial extends Simulation {

    protected Box box;
    protected IPotential2 potential1, potential2, potential12;
    protected IntegratorMC integrator;
    protected MCMoveAtom mcMoveAtom;
    protected MCMoveInsertDelete mcMoveInsertDelete;
    protected MCMoveGeometricCluster mcMoveGeometricCluster;
    protected SpeciesGeneral species1, species2;


    /**
     * @param numAtoms no. of solute atoms in the box
     * @param vf reservoir volume fraction of solvent
     * @param q size of solvent divided by size of solute
     * @param computeIdeal whether to compute histograms for ideal gas
     */
    public AshtonWildingOsmoticVirial(int numAtoms, double vf, double q, boolean computeIdeal, double L, double GCfreq){

        super(Space3D.getInstance());

        setRandom(new RandomMersenneTwister(1));

        species1 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species1);
        if (vf != 0) addSpecies(species2);

        double sigma1 = 1.0; //solute
        double sigma2 = q * sigma1; //solvent
        double sigma12 = (sigma1+sigma2)/2;

        box = new Box(new BoundaryRectangularPeriodic(space, L * sigma1), space);
        addBox(box);
        box.setNMolecules(species1,numAtoms);

        if(computeIdeal) {
            potential1 = null;
            potential2 = null;
            potential12 = null;
            System.out.println("P_ideal");
        }
        else{
            potential1 = P2HardSphere.makePotential(sigma1);
            potential2 = P2HardSphere.makePotential(sigma2);
            potential12 = P2HardSphere.makePotential(sigma12);
        }

        PotentialComputePair potentialMaster;
        NeighborManager neighborManager;
        if(!computeIdeal && vf != 0) {
            AtomType solvent = species2.getLeafType();
            AtomType solute = species1.getLeafType();
            neighborManager = new NeighborManagerCellMixed(getSpeciesManager(), box, 2, BondingInfo.noBonding(), solvent);
            potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
            potentialMaster.setPairPotential(solvent, solvent, potential2);
            potentialMaster.setPairPotential(solute, solvent, potential12);
        }
        else{
            neighborManager = new NeighborManagerSimple(box);
            potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        }
        potentialMaster.setPairPotential(species1.getLeafType(), species1.getLeafType(), potential1);

        integrator = new IntegratorMC(potentialMaster, random, 1.0, box);
        this.getController().addActivity(new ActivityIntegrate(integrator));
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        if(!computeIdeal && vf != 0) {
            mcMoveInsertDelete = new MCMoveInsertDelete(potentialMaster, random, space);
            mcMoveInsertDelete.setSpecies(species2);
            double mu = (8 * vf - 9 * vf * vf + 3 * vf * vf * vf) / Math.pow((1 - vf), 3) + Math.log(6 * vf / (Math.PI * Math.pow(sigma2, 3))); //Configurational chemical potential from Carnahan–Starling equation of state
            System.out.println("mu " + mu + " muig " + Math.log(6 * vf / (Math.PI * Math.pow(sigma2, 3))));
//              mu = Double.MAX_VALUE;
            mcMoveInsertDelete.setMu(mu);
            integrator.getMoveManager().addMCMove(mcMoveInsertDelete);
            mcMoveGeometricCluster = new MCMoveGeometricCluster(potentialMaster, neighborManager, space, random, integrator,
                    species1, potentialMaster.getPairPotentials());
            integrator.getMoveManager().addMCMove(mcMoveGeometricCluster);
        }

        System.out.println("vol: "+box.getBoundary().volume());

        if (vf != 0) {
            // add solvents to the boxes as though there are no solutes (and add solutes as though there
            // are no solvents).  then remove any solvents that overlap a solute
            int nSolvent = (int) (box.getBoundary().volume() * vf / (q * q * q * Math.PI / 6));
            box.setNMolecules(species2, nSolvent);
        }

//      integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box, species1);
        if (vf != 0) configuration.initializeCoordinates(box, species2);

        potentialMaster.init();
        if (!computeIdeal && vf != 0){
            final Set<IMolecule> overlaps = new HashSet<>();
            PotentialCallback overlapCheck = new PotentialCallback() {
                @Override
                public void pairCompute(int i, int j, Vector dr, double[] u012) {
                    if (u012[0] < Double.POSITIVE_INFINITY) return;
                    IAtom a = box.getLeafList().get(i);
                    if (a.getType().getSpecies() == species1) a = box.getLeafList().get(j);
                    overlaps.add(a.getParentGroup());
                }
            };

            potentialMaster.computeAll(false, overlapCheck);
            for (IMolecule m : overlaps) box.removeMolecule(m);
        }
        // GCfreq yields a simulation that spends roughly equal
        // CPU time doing GC moves and Displacement/Insert-Delete moves.  Set GCfreq to achieve desired result.

        if(!computeIdeal && vf != 0){
            integrator.getMoveManager().setFrequency(mcMoveGeometricCluster,2*GCfreq);
        }
    }

    public static void main(String[] args) throws IOException {
        simParams params = new simParams();

        if(args.length > 0){
            ParseArgs.doParseArgs(params, args);
        }
        else{
            params.numAtoms = 2;
            params.numSteps = 1000000;
            params.nBlocks = 100;
            params.vf = 0.05;
            params.computeIdeal = false;
            params.sizeRatio = 0.33333333333333333333;
            params.L = 3;
            params.GCfreq = 3;
            params.graphics = false;
        }
        int numAtoms = params.numAtoms;
        long numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double vf = params.vf;
        double q = params.sizeRatio;
        boolean computeIdeal = params.computeIdeal;
        double L = params.L;
        double GCfreq = params.GCfreq;
        boolean graphics = params.graphics;
        AccumulatorAverageFixed accNm = null;

        long numSamples = numSteps / numAtoms ;
        long samplesPerBlock = numSamples / nBlocks;
        if(samplesPerBlock == 0) samplesPerBlock = 1;

        System.out.println(numAtoms + " Atoms, "+ numSteps + " Steps" );
        System.out.println("Volume fraction: "+ vf);
        System.out.println("Size ratio: "+ q);
        System.out.println("nBlocks "+ nBlocks);
        System.out.println("system size: "+L+" sigma" );
        System.out.println("GC factor " + GCfreq);

        long t1 = System.currentTimeMillis();

        AshtonWildingOsmoticVirial sim = new AshtonWildingOsmoticVirial(numAtoms, vf, q, computeIdeal, L, GCfreq);

        if(graphics){
            final String appName = "Ashton-Wilding";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, appName, 3);

            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species1.getLeafType(), 1);
            if (vf != 0) ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            simGraphic.makeAndDisplayFrame(appName);
            return;
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps / 10));

        sim.integrator.getMoveManager().setEquilibrating(false);

        MeterRminSpecies meterRmin = new MeterRminSpecies(sim.space, sim.box, sim.species1);
        AccumulatorHistogram accRmin = new AccumulatorHistogram(new HistogramSimple(new DoubleRange(0, L)));
        DataPumpListener pumpRmin = new DataPumpListener(meterRmin,accRmin,numAtoms);
        sim.integrator.getEventManager().addListener(pumpRmin);

        if (vf != 0) {
            MeterNMolecules meterNMolecules = new MeterNMolecules();
            meterNMolecules.setSpecies(sim.species2);
            meterNMolecules.setBox(sim.box);
            accNm = new AccumulatorAverageFixed(samplesPerBlock);
            DataPumpListener pumpNm = new DataPumpListener(meterNMolecules, accNm);
            sim.integrator.getEventManager().addListener(pumpNm);
        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        double[] histRmin = accRmin.getHistograms().getHistogram();
        double[] r = accRmin.getHistograms().xValues();
        System.out.println("\nWriting probabilities to test.txt file\n");
        FileWriter writer = new FileWriter("test.txt");
        for (int i = 0; i < histRmin.length; i++) {
            writer.write(r[i]+" "+histRmin[i]+"\n");
        }
        writer.close();

        if (vf != 0) {
            double NmAvg = accNm.getData(AccumulatorAverage.AVERAGE).getValue(0);
            double NmErr = accNm.getData(AccumulatorAverage.ERROR).getValue(0);
            double NmCor = accNm.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);

            double VfBox = NmAvg * q * q * q * Math.PI / 6 / sim.box.getBoundary().volume();
            double VfErr = NmErr * q * q * q * Math.PI / 6 / sim.box.getBoundary().volume();

            System.out.print(String.format("No. of molecules avg: %13.6e  err: %11.4e   cor: % 4.2f\n", NmAvg, NmErr, NmCor));
            System.out.print(String.format("Volume frac avg: %13.6e  err: %11.4e\n", VfBox, VfErr));
        }

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+ (t2-t1)*0.001);

    }

    public static class simParams extends ParameterBase{
        public int numAtoms = 2;
        public long numSteps = 10000000;
        public int nBlocks = 100;
        public double vf = 0.0;
        public double sizeRatio = 0.33333333333333333333;
        public boolean computeIdeal = false;
        public double L = 3.0;
        public double GCfreq = 3;
        public boolean graphics = false;

    }
}
