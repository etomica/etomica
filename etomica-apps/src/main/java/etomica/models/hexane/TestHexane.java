/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

import etomica.action.PDBWriter;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorListenerGroupSeries;
import etomica.normalmode.*;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.virial.MCMoveClusterWiggleMulti;
/**
 * @author nancycribbin
 *  
 */

/*
 * We use a PotentialMaster, rather than a PotentialMasterNbr, so that we do not
 * need to deal with cells, which BoundaryDeformablePeriodic cannot deal with at
 * this time.
 * 
 */

public class TestHexane extends Simulation {

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;

    public Box box;

    public BoundaryDeformablePeriodic bdry;
    public BravaisLattice lattice;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    
    public MCMoveClusterWiggleMulti crank; 
//    public MCMoveReptate snake;
    public MCMoveMolecule moveMolecule;
    public CBMCGrowSolidHexane growMolecule;
    public MCMoveRotateMolecule3D rot;
    public MCMoveMoleculeCoupled coupledMove;
    public MCMoveCombinedCbmcTranslation cctMove;
    
//    public PairIndexerMolecule pri;

    
    public TestHexane(Space _space, double dens, int xCells, int yCells, int zCells) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(_space);

        PotentialMaster potentialMaster = new PotentialMaster();
        int chainLength = 6;
        //One molecule per cell
        int numAtoms = xCells * yCells * zCells * chainLength;
        primitive = new PrimitiveHexane(space);
        // close packed density is 0.4165783882178116
        // Monson reports data for 0.373773507616 and 0.389566754417
        primitive.scaleSize(Math.pow(0.4165783882178116/dens,1.0/3.0));
        lattice = new BravaisLattice(primitive);

        //This is the factor that multiples by the range of the potential in
        // order to define the area/volume in which neighbors are searched for.
        //This becomes the bond delta, which is the percentage the bond can
        // stretch, and I assume compress.
        double neighborRangeFac = 1.2;

        double bondFactor = 0.4;

        SpeciesHexane species = new SpeciesHexane(space);
        addSpecies(species);
        int[] nCells = new int[]{xCells, yCells, zCells};
        bdry = new BoundaryDeformableLattice(primitive, nCells);
        box = new Box(bdry, space);
        addBox(box);
        box.setNMolecules(species, xCells * yCells * zCells);
//        config.initializeCoordinates(box);
        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0, box);
        
        moveMolecule = new MCMoveMolecule(potentialMaster, getRandom(), space,
                0.1, 1);
        // 0.025 for translate, 0.042 for rotate for rho=0.3737735
        moveMolecule.setStepSize(0.024);        
        integrator.getMoveManager().addMCMove(moveMolecule);
        ((MCMoveStepTracker)moveMolecule.getTracker()).setNoisyAdjustment(true);
        
//        moveVolume = new MCMoveVolume(this, potentialMaster);
////        moveVolume = new MCMoveVolume(potentialMaster, getRandom(), pressure);
//        moveVolume.setBox(box);
//        integrator.getMoveManager().addMCMove(moveVolume);
        
        crank = new MCMoveClusterWiggleMulti(potentialMaster, getRandom(), 0.20, 6, space);
    
//        snake = new MCMoveReptate(potentialMaster, getRandom(), 0.4, 3.0, true);
//        snake.setBox(box);
//        integrator.getMoveManager().addMCMove(snake);
        
        rot = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
        rot.setBox(box);
        rot.setStepSize(0.042);
        integrator.getMoveManager().addMCMove(rot);
        ((MCMoveStepTracker)rot.getTracker()).setNoisyAdjustment(true);
        
        growMolecule = new CBMCGrowSolidHexane(potentialMaster,
                getRandom(), space, integrator, box, species, 20);
        growMolecule.setBox(box);
        integrator.getMoveManager().addMCMove(growMolecule);

        coupledMove = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
        integrator.getMoveManager().addMCMove(coupledMove);
        
        cctMove = new MCMoveCombinedCbmcTranslation(potentialMaster, growMolecule, getRandom(), space);
        cctMove.setBox(box);
        integrator.getMoveManager().addMCMove(cctMove);
        
        // nan we're going to need some stuff in there to set the step sizes and
        // other stuff like that.

        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(2000000);
        getController().addAction(activityIntegrate);
            
        //nan The box size we want is 5.72906360610622 by 11.21417818673970 by
        // 7.30591061708510
        //nan this is where the squared, unsquared box stuff comes in.
        //makes the density 0.41657 per Dr. Monson's comment in e-mail.
//        defaults.boxSize = 7.018;
//        defaults.boxSize = 100;

        //INTERMOLECULAR POTENTIAL STUFF

        //This potential is the intermolecular potential between atoms on
        // different molecules. We use the class "Potential" because we are
        // reusing the instance as we define each potential.
        Potential potential = new P2HardSphere(space);
        
        //here, we add the species to the PotentialMaster, using types.
        //The PotentialMaster generates a group potential and automatically
        // does a lot of the stuff which we have to do for the intramolecular
        // potential manually.
        AtomType sphereType = species.getLeafType();

        //Add the Potential to the PotentialMaster
        potentialMaster.addPotential(potential, new AtomType[]{sphereType,
                sphereType });
        
        coupledMove.setPotential(potentialMaster.getPotential(new ISpecies[] {
                species, species }  ));

        //Initialize the positions of the atoms.
        coordinateDefinition = new CoordinateDefinitionHexane(this, box, primitive, species, space);
        coordinateDefinition.initializeCoordinates(nCells);
//        WriteConfiguration writer = new WriteConfiguration();
//        writer.setBox(box);
//        writer.setDoApplyPBC(false);
//        writer.setConfName("hexanePure");
//        writer.actionPerformed();
        
        integrator.setBox(box);
       
    }

    public static void main(String[] args) {
        int xLng = 4;
        int yLng = 4;
        int zLng = 3;
        long nSteps = 1000000;
        // Monson reports data for 0.373773507616 and 0.389566754417
        double density = 0.349942899;
        double den = 0.35;
        boolean graphic = true;
  
        //spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        TestHexane sim = new TestHexane(Space3D.getInstance(), density, xLng, yLng, zLng);

        System.out.println("Happy Goodness!!");

        if (graphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.makeAndDisplayFrame();
        } else {
            long time = System.currentTimeMillis();
            System.out.println(time);
            long time1;
            long time2;
            
            //parse arguments
            //filename is element 0
            String filename = "nm_hex_";
                
            if(args.length > 0){
                filename = args[0];
            }
            if(args.length > 1){
                nSteps = Long.parseLong(args[1]);
            }
            if(args.length > 2){
                den = Double.parseDouble(args[2]);
                if(den == 0.37) {density = 0.373773507616;}
                if(den == 0.40) {density = 0.389566754417;}
                if(den == 0.35) {density = 0.349942899;}
            }
            if(args.length > 5){
                xLng = Integer.parseInt(args[3]);
                yLng = Integer.parseInt(args[4]);
                zLng = Integer.parseInt(args[5]);
            }
            filename = filename + nSteps + "_" + (int)(den*100) + "_" + xLng 
                + "_" + yLng + "_" + zLng;
            
            System.out.println(filename);
            System.out.println("Hexane simulation");
            System.out.println("Number of steps = " + nSteps);
            System.out.println("Density = " + density);
            System.out.println("Number of cells/molecules in the X direction = " + xLng);
            System.out.println("Number of cells/molecules in the Y direction = " + yLng);
            System.out.println("Number of cells/molecules in the Z direction = " + zLng);
            System.out.println("Total number of molecules = " + xLng*yLng*zLng);
            double volume = sim.bdry.volume();
            System.out.println("volume =  "+ volume);
            
            PrimitiveHexane primitive = (PrimitiveHexane)sim.lattice.getPrimitive();
            // primitive doesn't need scaling.  The boundary was designed to be commensurate with the primitive
            WaveVectorFactorySimple waveVectorFactory = new WaveVectorFactorySimple(primitive, sim.space);
            // we need to set this up now even though we don't use it during equilibration so that
            // the meter can grab the lattice points
            MeterNormalMode meterNormalMode = new MeterNormalMode();
            meterNormalMode.setWaveVectorFactory(waveVectorFactory);
            meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
            meterNormalMode.setBox(sim.box);

//            BoxInflateDeformable pid = new BoxInflateDeformable(sim.getSpace());
//            MeterPressureByVolumeChange meterPressure = new MeterPressureByVolumeChange(sim.getSpace(), pid);
//            meterPressure.setIntegrator(sim.integrator);
//            meterPressure.setX(-0.001, 0.001, 20);
//            AccumulatorAverageFixed pressureAccumulator = new AccumulatorAverageFixed();
//            DataPump pressureManager = new DataPump(meterPressure, pressureAccumulator);
//            pressureAccumulator.setBlockSize(50);
//            sim.integrator.addIntervalAction(pressureManager);
//            sim.integrator.setActionInterval(pressureManager, 10);
         
//            sim.activityIntegrate.setMaxSteps(nSteps/10);
            sim.activityIntegrate.setMaxSteps(30000);
            sim.getController().actionPerformed();
            System.out.println("equilibration finished: " + sim.activityIntegrate.getMaxSteps() +"  steps");
            sim.getController().reset();
            
//            IAction progressReport = new IAction() {
//                public void actionPerformed() {
//                    System.out.print(sim.integrator.getStepCount()+" steps: ");
//                }
//            };
            
            
            ((MCMoveStepTracker)sim.moveMolecule.getTracker()).setTunable(false);
            ((MCMoveStepTracker)sim.rot.getTracker()).setTunable(false);
           
            IntegratorListenerGroupSeries listenerGroup = new IntegratorListenerGroupSeries();
            
            IntegratorListenerAction meterListener = new IntegratorListenerAction(meterNormalMode);
            meterListener.setInterval(10000);
            listenerGroup.addListener(meterListener);
            
            DataGroup normalModeData = (DataGroup)meterNormalMode.getData();
//            normalModeData.TE(1.0/(sim.box.getSpeciesMaster().moleculeCount()*meterNormalMode.getCallCount()));
            normalModeData.TE(1.0/sim.box.getMoleculeList().getMoleculeCount());
            int normalDim = meterNormalMode.getCoordinateDefinition().getCoordinateDim();
            
            Vector[] waveVectors = waveVectorFactory.getWaveVectors();
            double[] coefficients = waveVectorFactory.getCoefficients();
            
            WriteS sWriter = new WriteS(sim.space);
            sWriter.setFilename(filename);
            sWriter.setOverwrite(true);
            sWriter.setMeter(meterNormalMode);
            sWriter.setWaveVectorFactory(waveVectorFactory);
            sWriter.setTemperature(sim.integrator.getTemperature());
            
            IntegratorListenerAction sWriterListener = new IntegratorListenerAction(sWriter);
            sWriterListener.setInterval((int)nSteps/10);
            listenerGroup.addListener(sWriterListener);
            sim.integrator.getEventManager().addListener(listenerGroup);
            
            sim.activityIntegrate.setMaxSteps(nSteps);
            sim.getController().actionPerformed();

            //Calls the harmonic calculator to be sure that no data is lost
            sWriter.actionPerformed();
            
            //Write out the final configurations for further use.
            PDBWriter pdbWriter = new PDBWriter();
            pdbWriter.setBox(sim.box);
            pdbWriter.setFileName("calcHex.pdb");
            pdbWriter.actionPerformed();

            WriteConfiguration writer = new WriteConfiguration(sim.space);
            writer.setBox(sim.box);
            writer.setDoApplyPBC(false);
            writer.setConfName("hexane");
            writer.actionPerformed();
            time1 = System.currentTimeMillis();   
            
//            double avgPressure = 0.0;  
//            int leng = 20;
//            double[] pressies = new double[leng];
//            double[] lnXs = new double[leng];
//            double[] scalingFactors = new double[leng];
//            double[] volumes = new double[leng];
//
//            lnXs = ((DataDoubleArray)((DataGroup)pressureAccumulator.getData()).getData(StatType.AVERAGE.index)).getData();
//            
//            for(int i = 0; i < leng; i++){
//                scalingFactors[i] = ((DataDoubleArray)meterPressure.getScalingDataSource().getData()).getValue(i);
//                lnXs[i] = Math.log(lnXs[i]);
//                volumes[i] = volume*scalingFactors[i];
//                pressies[i] = lnXs[i]/volumes[i];
//            }
//            
//            System.out.println("volume =  "+ volume);
//            System.out.println("lnXs");
//            for(int i = 0; i < leng; i++){
//                System.out.println(lnXs[i]);
//            }
//            System.out.println("volumes");
//            for(int i = 0; i < leng; i++){
//                System.out.println(volumes[i]);
//            } 
//            System.out.println("scaling factors");
//            for (int i = 0; i < leng; i++){
//                System.out.println(scalingFactors[i]);
//            }
//            
//            avgPressure = ((DataDoubleArray)((DataGroup)pressureAccumulator.getData()).getData(StatType.AVERAGE.index)).getValue(0);
//            System.out.println("Avg Pres = "+ avgPressure);
//            time2 = System.currentTimeMillis();
//            
//            System.out.println("start  " + time);
//            System.out.println("simulation  " + time1);
//            System.out.println("data collection  " + time2);
            
//            System.out.println(sim.integrator.meterPE.getDataAsScalar());
            
            System.out.println("Fini.");
        }
    }
}
