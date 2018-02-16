package etomica.osmoticvirial;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
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
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2HardSphere;
import etomica.potential.P2Ideal;
import etomica.potential.Potential2;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

/**
 * Implements Ashton and Wilding method for calculation of osmotic virial coefficient (B2, B3)
 * for Hard-Sphere model as described in the paper DOI: 10.1063/1.4883718.
 */
public class AshtonWildingOsmoticVirial extends Simulation {

    protected Box box;
    protected Potential2 potential1, potential2, potential12;
    protected Controller controller;
    protected IntegratorMC integrator;
    protected MCMoveAtom mcMoveAtom;
    protected MCMoveInsertDelete mcMoveInsertDelete;
    protected MCMoveGeometricCluster mcMoveGeometricCluster;
    protected SpeciesSpheresMono species1, species2;
    protected ActivityIntegrate activityIntegrate;

    /**
     * @param numAtoms no. of solute atoms in the box
     * @param vf reservoir volume fraction of solvent
     * @param q size of solvent divided by size of solute
     * @param computeIdeal whether to compute histograms for ideal gas
     */
    public AshtonWildingOsmoticVirial(int numAtoms, double vf, double q, boolean computeIdeal){

        super(Space3D.getInstance());
        setRandom(new RandomMersenneTwister(1));
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, space);

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
        System.out.println("vol: "+box.getBoundary().volume());

        mcMoveInsertDelete.setSpecies(species2);
        double mu = (8*vf-9*vf*vf+3*vf*vf*vf) / Math.pow((1-vf),3)+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))); //Configurational chemical potential from Carnahan–Starling equation of state
        System.out.println("mu "+ mu+" muig "+Math.log(6*vf/(Math.PI*Math.pow(sigma2,3))));
       // mu = Double.MAX_VALUE;
        mcMoveInsertDelete.setMu(mu);

        box.setNMolecules(species1, numAtoms);

        if(computeIdeal) {
            potential1 = new P2Ideal(space);
            potential2 = new P2Ideal(space);
            potential12 = new P2Ideal(space);
            ((P2Ideal) potential1).setRange(1);
            System.out.println("P_ideal");
        }
        else{
            potential1 = new P2HardSphere(space, sigma1, false);
            potential2 = new P2HardSphere(space, sigma2, false);
//            potential2 = new P2Ideal(space);
            potential12 = new P2HardSphere(space, sigma12, false);
//            System.out.println("AO");
            mcMoveGeometricCluster = new MCMoveGeometricCluster(potentialMaster, space, random, integrator, species1);
            integrator.getMoveManager().addMCMove(mcMoveGeometricCluster);
        }

        potentialMaster.setCellRange(3);
        potentialMaster.setRange(potential1.getRange());
        potentialMaster.setPotentialHard(true);

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

        if(args.length > 0){
            ParseArgs.doParseArgs(params, args);
        }
        else{
            params.numAtoms = 3;
            params.numSteps = 100000;
            params.nBlocks = 100;
            params.vf = 0.1;
            params.computeIdeal = false;
            params.sizeRatio = 0.2;
        }
        int numAtoms = params.numAtoms;
        long numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double vf = params.vf;
        double q = params.sizeRatio;
        boolean computeIdeal = params.computeIdeal;
        boolean graphics = false;
        AccumulatorAverageFixed accNm = null;

        long numSamples = numSteps / numAtoms ;
        long samplesPerBlock = numSamples / nBlocks;
        if(samplesPerBlock == 0) samplesPerBlock = 1;

        System.out.println(numAtoms + " Atoms, "+ numSteps + " Steps" );
        System.out.println("Volume fraction: "+ vf);
        System.out.println("Size ratio: "+ q);
        System.out.println("nBlocks "+ nBlocks);

        long t1 = System.currentTimeMillis();

        AshtonWildingOsmoticVirial sim = new AshtonWildingOsmoticVirial(numAtoms, vf, q, computeIdeal);

        if(graphics){
            final String appName = "Ashton-Wilding";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, appName, 3, sim.getSpace(), sim.getController());

            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species1.getLeafType(), 1);
            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            simGraphic.makeAndDisplayFrame(appName);
            return;
        }
        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);

        MeterRminSpecies meterRmin = new MeterRminSpecies(sim.space, sim.box, sim.species1);
        AccumulatorHistogram accRmin = new AccumulatorHistogram(new HistogramSimple(new DoubleRange(0, 4*sim.potential1.getRange())));
        DataPumpListener pumpRmin = new DataPumpListener(meterRmin,accRmin,numAtoms);
        sim.integrator.getEventManager().addListener(pumpRmin);

        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.species2);
        meterNMolecules.setBox(sim.box);
        accNm = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpNm = new DataPumpListener(meterNMolecules, accNm);
        sim.integrator.getEventManager().addListener(pumpNm);

        sim.getController().actionPerformed();
        double[] histRmin = accRmin.getHistograms().getHistogram();
        double[] r = accRmin.getHistograms().xValues();

        for (int i = 0; i < histRmin.length; i++) {
            System.out.println(r[i]+" "+histRmin[i]);
        }

        double NmAvg = accNm.getData(AccumulatorAverage.AVERAGE).getValue(0);
        double NmErr = accNm.getData(AccumulatorAverage.ERROR).getValue(0);
        double NmCor = accNm.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);

        double VfBox = NmAvg*q*q*q*Math.PI/6/sim.box.getBoundary().volume();
        double VfErr = NmErr*q*q*q*Math.PI/6/sim.box.getBoundary().volume();

        System.out.print(String.format("No. of molecules avg: %13.6e  err: %11.4e   cor: % 4.2f\n", NmAvg, NmErr, NmCor));
        System.out.print(String.format("Volume frac avg: %13.6e  err: %11.4e\n", VfBox, VfErr));

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+ (t2-t1)*0.001);

    }

    public static class simParams extends ParameterBase{
        public int numAtoms = 20;
        public long numSteps = 10000;
        public int nBlocks = 100;
        public double vf = 0.2;
        public double sizeRatio = 0.2;
        public boolean computeIdeal = false;

    }
}
