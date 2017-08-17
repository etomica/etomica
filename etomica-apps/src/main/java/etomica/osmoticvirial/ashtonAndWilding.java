package etomica.osmoticvirial;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.histogram.HistogramSimple;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.DoubleRange;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Implements Ashton and Wilding method for calculation of gas virial coefficient for Hard-Sphere model
 * as described in the paper.
 * Created by aksharag on 07-11-2017.
 */


public class ashtonAndWilding extends Simulation {

    public Box box;
    public Potential2 potential1;
    public Controller controller;
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species1;
    public ActivityIntegrate activityIntegrate;


    public ashtonAndWilding(int numAtoms, double density, boolean computeIdeal, boolean computeDepletion){

        super(Space3D.getInstance());
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, space);

        integrator = new IntegratorMC(this, potentialMaster);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);

        integrator.getMoveManager().addMCMove(mcMoveAtom);

        double sigma1 = 1.0;

        box = new Box(space);
        addBox(box);
        species1 = new SpeciesSpheresMono(this, space);
        addSpecies(species1);
        box.setNMolecules(species1, numAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        if(!computeDepletion) potential1 = new P2HardSphere(space, sigma1, false);
        if(computeDepletion) {
            potential1 = new PotentialDepletion(space, "/usr/users/aksharag/osmoticvirials/ashtonWilding/additiveHS/mi+me/more steps/10e10/g3");
        }

        potentialMaster.setCellRange(3);
        potentialMaster.setRange(potential1.getRange());

        AtomType leafType1 = species1.getLeafType();

       if(!computeIdeal) potentialMaster.addPotential(potential1, new AtomType[]{leafType1, leafType1});
       if(computeIdeal) System.out.println("P_ideal");

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
            params.numSteps = 1000000;
            params.nBlocks = 100;
            params.density = 0.01;
            params.computeIdeal = false;
        }

        int numAtoms = params.numAtoms;
        int numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double density = params.density;
        boolean computeIdeal = params.computeIdeal;
        boolean graphics = false;
        boolean computeDep = true;

        long numSamples = numSteps / numAtoms;
        long samplesPerBlock = numSamples / nBlocks;
        if(samplesPerBlock == 0) samplesPerBlock = 1;

        System.out.println(numAtoms + " atoms, "+ numSteps + " steps" );
        System.out.println("density: "+ density);
        System.out.println("nBlocks "+ nBlocks);

        long t1 = System.currentTimeMillis();

        ashtonAndWilding sim = new ashtonAndWilding(numAtoms, density, computeIdeal, computeDep);

        if(graphics){
            final String appName = "Ashton-Wilding";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, appName, 3, sim.getSpace(), sim.getController());

            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species1.getLeafType(), 1);
            simGraphic.makeAndDisplayFrame(appName);

            return;
        }

        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);

        MeterRmin meterRmin = new MeterRmin(sim.space, sim.box);
        AccumulatorHistogram accRmin = new AccumulatorHistogram(new HistogramSimple(new DoubleRange(0, 4*sim.potential1.getRange())));
        DataPumpListener pumpRmin = new DataPumpListener(meterRmin,accRmin,numAtoms);
        sim.integrator.getEventManager().addListener(pumpRmin);
        sim.getController().actionPerformed();
        double[] histRmin = accRmin.getHistograms().getHistogram();
        double[] r = accRmin.getHistograms().xValues();

        for (int i = 0; i < histRmin.length; i++)

        {
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
        public boolean computeIdeal = true;

    }
}
