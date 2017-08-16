package etomica.osmoticvirial;

import etomica.action.ActionIntegrate;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterWidomInsertion;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Created by aksharag on 6/16/17.
 */
public class fepLJ extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
  //  public MCMoveInsertDelete mcMoveInsertDelete : for converting nvt to mu-v-t. ; we are adding and removing solvent molecules (changing N)
    public SpeciesSpheresMono species1; //solvent
    public SpeciesSpheresMono species2; //solute
    public Box box;
    public P2LennardJones potential1, potential2, potential12;
    public Controller controller;
    public ActivityIntegrate activityIntegrate;

    public fepLJ(int numAtoms, int numSteps, double temp, double density, double sigma2, double epsilon2, boolean computez2){
        super(Space3D.getInstance());
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this,space);

        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setTemperature(temp);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
      //  mcMoveInsertDelete = new MCMoveInsertDelete(potentialMaster, random, space);

        integrator.getMoveManager().addMCMove(mcMoveAtom);
      //  integrator.getMoveManager().addMCMove(mcMoveInsertDelete);

        double sigma1 = 1.0;
        double epsilon1 = 1.0;
        double sigma12 = (sigma1+sigma2)/2;
        double epsilon12 = Math.sqrt(epsilon1*epsilon2) ;

        species1 = new SpeciesSpheresMono(this, space);
        species2 = new SpeciesSpheresMono(this, space);
        addSpecies(species1);
        addSpecies(species2);
        box = new Box(space);
        addBox(box);

       // mcMoveInsertDelete.setSpecies(species1);
        //mcMoveInsertDelete.setMu(); //TODO


        box.setNMolecules(species1,numAtoms);

        BoxInflate inflater = new BoxInflate(box,space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        if (computez2){box.setNMolecules(species2,1);}


        potential1 = new P2LennardJones(space, sigma1, epsilon1);
        potential2 = new P2LennardJones(space, sigma2, epsilon2);
        potential12 = new P2LennardJones(space, sigma12, epsilon12);

        double truncationRadius1 = 3.0*sigma1;
        double truncationRadius2 = 3.0*sigma1;
        double truncationRadius12 = 3.0*sigma1;

        if(truncationRadius1>0.5*box.getBoundary().getBoxSize().getX(0)){
            throw new RuntimeException(" Truncation radius is too large. Max allowed is:"+ 0.5*box.getBoundary().getBoxSize().getX(0));
        }

        if(truncationRadius2>0.5*box.getBoundary().getBoxSize().getX(0)){
            throw new RuntimeException(" Truncation radius is too large. Max allowed is:"+ 0.5*box.getBoundary().getBoxSize().getX(0));
        }

        potentialMaster.setCellRange(3);

        P2SoftSphericalTruncated potentialTruncated11 = new P2SoftSphericalTruncated(space, potential1, truncationRadius1);
        potentialMaster.setRange(potentialTruncated11.getRange());
        AtomType leafType1 = species1.getLeafType();
        potentialMaster.addPotential(potentialTruncated11, new AtomType[]{leafType1, leafType1});

        P2SoftSphericalTruncated potentialTruncated12 = new P2SoftSphericalTruncated(space, potential12, truncationRadius12);
        potentialMaster.setRange(potentialTruncated12.getRange());
        AtomType leafType2 = species2.getLeafType();
        potentialMaster.addPotential(potentialTruncated12, new AtomType[]{leafType1, leafType2});

        P2SoftSphericalTruncated potentialTruncated22 = new P2SoftSphericalTruncated(space, potential2, truncationRadius2);
        potentialMaster.setRange(potentialTruncated22.getRange());
        potentialMaster.addPotential(potentialTruncated22, new AtomType[]{leafType2, leafType2});

        // potentialMaster.setRange should pass the max of the three truncation radii as the parameter


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
            params.numAtoms = 500;
            params.numSteps = 500000;
            params.nBlocks = 1000;
            params.temp = 2;
            params.density = 0.6;
            params.sigma2 = 1.5;
            params.eps2 = 1;
            params.computez2 = true;
        }

        int numAtoms = params.numAtoms;
        int numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double temp = params.temp;
        double density = params.density;
        double sigma2 = params.sigma2;
        double eps2 = params.eps2;
        boolean computez2 = params.computez2;
        boolean graphics = true;

        long numSamples = numSteps/numAtoms;
        long samplesPerBlock = numSamples/nBlocks;
        if (samplesPerBlock == 0) samplesPerBlock = 1;

        if(computez2){
            System.out.println("**z2_z1**");
        }
        else{
            System.out.println("**z1_z0**");
        }

        System.out.println(numAtoms+" atoms, "+numSteps+" steps");
        System.out.println("density: "+density);
        System.out.println("temperature: "+temp);
        System.out.println("sigma2: "+sigma2);
        System.out.println(nBlocks+" blocks");

        long t1 = System.currentTimeMillis();


        fepLJ sim = new fepLJ(numAtoms, numSteps, temp, density, sigma2, eps2, computez2);

        System.out.println("box length "+sim.box.getBoundary().getBoxSize());

        if (graphics) {
            final String APP_NAME = "SimLJ";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());

            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species2.getLeafType(), sigma2);
            simGraphic.makeAndDisplayFrame(APP_NAME);

            MeterWidomInsertion meterinsert = new MeterWidomInsertion(sim.space,sim.getRandom());
            //meterinsert.setNInsert(50);
            meterinsert.setSpecies(sim.species2);
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
        meterinsert.setSpecies(sim.species2);
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
        public double temp = 2;
        public double density = 0.2;
        public double sigma2 = 2.0;
        public double eps2 = 1.0;
        public boolean computez2 = false;
    }
}
