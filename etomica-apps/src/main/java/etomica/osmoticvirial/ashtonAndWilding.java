package etomica.osmoticvirial;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomsetIteratorBoxDependent;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Created by Lenovo on 07-11-2017.
 */


public class ashtonAndWilding extends Simulation {

    public Box box;
    public P2HardSphere potential1;
    public Controller controller;
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species1;
    public ActivityIntegrate activityIntegrate;
    public AtomsetIteratorBoxDependent iterator;


    public ashtonAndWilding(double temp, int numAtoms, double density){

        super(Space3D.getInstance());
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, space);

        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setTemperature(temp);
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

        potential1 = new P2HardSphere(space, sigma1, false);

        potentialMaster.setCellRange(3);
        potentialMaster.setRange(potential1.getRange());

        AtomType leafType1 = species1.getLeafType();

        potentialMaster.addPotential(potential1, new AtomType[]{leafType1, leafType1});

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);

        integrator.setBox(box);
        potentialMaster.getNbrCellManager(box).assignCellAll();

        iterator = new ApiLeafAtoms();

    }

    public static void main(String[] args){
        simParams params = new simParams();

        if(args.length > 0){
            ParseArgs.doParseArgs(params, args);
        }
        else{
            params.numAtoms = 500;
            params.numSteps = 50000;
            params.nBlocks = 100;
            params.density = 0.5;
            params.temp = 2.0;
        }

        int numAtoms = params.numAtoms;
        int numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double density = params.density;
        double temp = params.temp;
        boolean graphics = true;

        long numSamples = numSteps / numAtoms;
        long samplesPerBlock = numSamples / nBlocks;
        if(samplesPerBlock == 0) samplesPerBlock = 1;

        System.out.println(numAtoms + " atoms, "+ numSteps + " steps" );
        System.out.println("density: "+ density);
        System.out.println("temperature "+temp);
        System.out.println(" nBlocks"+ nBlocks);

        long t1 = System.currentTimeMillis();

        ashtonAndWilding sim = new ashtonAndWilding(temp, numAtoms, density);

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

        double min = 100;
        final Vector dr = sim.space.makeVector();
        AtomType type1, type2;

        for (IAtomList pair = sim.iterator.next(); pair != null;
             pair = sim.iterator.next()) {
            if (type1 != null && (pair.getAtom(0).getType() != type1 || pair.getAtom(1).getType() != type2)) continue;
           // dr.Ev1Mv2(pair.getAtom(1).getPosition(),pair.getAtom(0).getPosition());
           // boundary.nearestImage(dr);
            double r2 = dr.squared();
            if(min > r2) min = r2;
        }



        long t2 = System.currentTimeMillis();
        System.out.println("time: "+ (t2-t1)*0.001);

    }

    public static class simParams extends ParameterBase{
        public int numAtoms = 500;
        public int numSteps = 10000;
        public int nBlocks = 100;
        public double temp = 1.0;
        public double density = 0.1;

    }
}
