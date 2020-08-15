package etomica.osmoticvirial;

import etomica.action.BoxInflate;

import etomica.action.activity.ActivityIntegrate2;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPumpListener;
import etomica.data.histogram.HistogramSimple;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.DoubleRange;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Extends Ashton and Wilding approach for calculation of gas virial coefficient
 * to Lennard-Jones model.
 */
public class AshtonWildingLJ extends Simulation {

    protected Box box;
    protected P2LennardJones potential;
    protected IntegratorMC integrator;
    protected MCMoveAtom mcMoveAtom;
    protected SpeciesSpheresMono species1;
    

    /**
     * @param numAtoms number of atoms in the box
     * @param density density of the box
     * @param computeIdeal whether to compute histograms for ideal gas
     */
    public AshtonWildingLJ(int numAtoms, double density, boolean computeIdeal){
        super(Space3D.getInstance());

        double sigma1 = 1.0;
        box = new Box(space);
        addBox(box);
        species1 = new SpeciesSpheresMono(this, space);
        addSpecies(species1);
        box.setNMolecules(species1, numAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, space);
        integrator = new IntegratorMC(this, potentialMaster, box);
        this.getController2().addActivity(new ActivityIntegrate2(integrator));
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        potential = new P2LennardJones(space, sigma1, 1);
        double truncationRadius1 = 4.0*sigma1;
        if(truncationRadius1>0.5*box.getBoundary().getBoxSize().getX(0)){
            throw new RuntimeException(" Truncation radius is too large. Max allowed is:"+ 0.5*box.getBoundary().getBoxSize().getX(0));
        }
        potentialMaster.setCellRange(3);
        P2SoftSphericalTruncated potentialTruncated11 = new P2SoftSphericalTruncated(space, potential, truncationRadius1);
        potentialMaster.setRange(potentialTruncated11.getRange());
        AtomType leafType1 = species1.getLeafType();

        if(!computeIdeal) potentialMaster.addPotential(potentialTruncated11, new AtomType[]{leafType1, leafType1});
        if(computeIdeal) System.out.println("P_ideal");

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
        potentialMaster.getNbrCellManager(box).assignCellAll();
    }

    public static void main(String[] args){
        simParams params = new simParams();

        if(args.length > 0){
            ParseArgs.doParseArgs(params, args);
        }
        else{
            params.numAtoms = 2;
            params.numSteps = 1000000;
            params.nBlocks = 100;
            params.density = 0.003;
            params.computeIdeal = true;
            params.temp = 1;
        }
        int numAtoms = params.numAtoms;
        int numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double density = params.density;
        double temp = params.temp;
        boolean computeIdeal = params.computeIdeal;
        boolean graphics = false;

        long numSamples = numSteps / numAtoms;
        long samplesPerBlock = numSamples / nBlocks;
        if(samplesPerBlock == 0) samplesPerBlock = 1;

        System.out.println(numAtoms + " atoms, "+ numSteps + " steps" );
        System.out.println("density: "+ density);
        System.out.println("nBlocks "+ nBlocks);

        long t1 = System.currentTimeMillis();

        AshtonWildingLJ sim = new AshtonWildingLJ(numAtoms, density, computeIdeal);

        if(graphics){
            final String appName = "Ashton-Wilding";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, appName, 3);

            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species1.getLeafType(), 1);
            simGraphic.makeAndDisplayFrame(appName);
            return;
        }

        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps/10);
sim.getController().reset();
sim.integrator.getMoveManager().setEquilibrating(false);
        sim.integrator.setTemperature(temp);

        MeterRmin meterRmin = new MeterRmin(sim.space, sim.box);
        AccumulatorHistogram accRmin = new AccumulatorHistogram(new HistogramSimple(new DoubleRange(0, 0.9*sim.box.getBoundary().getBoxSize().getX(0))));
        DataPumpListener pumpRmin = new DataPumpListener(meterRmin,accRmin,numAtoms);
        sim.integrator.getEventManager().addListener(pumpRmin);
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps);
        double[] histRmin = accRmin.getHistograms().getHistogram();
        double[] r = accRmin.getHistograms().xValues();

        for (int i = 0; i < histRmin.length; i++) {
            System.out.println(r[i]+" "+histRmin[i]);
        }
        long t2 = System.currentTimeMillis();
        System.out.println("time: "+ (t2-t1)*0.001);
    }

    public static class simParams extends ParameterBase{
        public int numAtoms = 20;
        public int numSteps = 10000;
        public int nBlocks = 100;
        public double density = 0.1;
        public double temp = 1;
        public boolean computeIdeal = true;
    }
}
