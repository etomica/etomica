/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.*;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Pixel;

import java.util.ArrayList;

public class SimHarmonicHexane extends Simulation {

    private static final String APP_NAME = "Sim Harmonic";
    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;
    public Box box;
    public BoundaryDeformablePeriodic bdry;
    public BravaisLattice lattice;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    public NormalModes normalModes;
    public double fudge;

    
    public SimHarmonicHexane(Space _space, double dens, int xCells, int yCells, 
            int zCells, String filename, double harmonicFudge) {
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster();
        int chainLength = 6;
        //One molecule per cell
        int numAtoms = xCells * yCells * zCells * chainLength;
        fudge = harmonicFudge;

        SpeciesHexane species = new SpeciesHexane(space);
        addSpecies(species);

        primitive = new PrimitiveHexane(space);
        // close packed density is 0.4165783882178116
        // Monson reports data for 0.373773507616 and 0.389566754417
        primitive.scaleSize(Math.pow(0.4165783882178116 / dens, 1.0 / 3.0));
        lattice = new BravaisLattice(primitive);

        int[] nCells = new int[]{xCells, yCells, zCells};
        bdry = new BoundaryDeformableLattice(primitive, nCells);
        box = this.makeBox(bdry);
        box.setNMolecules(species, xCells * yCells * zCells);
//        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0);

        integrator = new IntegratorMC(this, null, box);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        MCMoveHarmonic harm = new MCMoveHarmonic(getRandom());
        integrator.getMoveManager().addMCMove(harm);

        coordinateDefinition = new CoordinateDefinitionHexane(this, box, primitive,
                species, space);
        coordinateDefinition.initializeCoordinates(nCells);

        normalModes = new NormalModesFromFile(filename, 3);
        normalModes.setHarmonicFudge(harmonicFudge);

        WaveVectorFactory waveFactory = normalModes.getWaveVectorFactory();
        waveFactory.makeWaveVectors(box);

        harm.setOmegaSquared(normalModes.getOmegaSquared());
        harm.setEigenVectors(normalModes.getEigenvectors());
        harm.setWaveVectors(waveFactory.getWaveVectors());
        harm.setWaveVectorCoefficients(waveFactory.getCoefficients());
        harm.setCoordinateDefinition(coordinateDefinition);
        harm.setBox(box);
    }

    public static void main(String[] args) {
//        int numMolecules = 144; //144

        //set up default values
        int xLng = 4;
        int yLng = 6;
        int zLng = 6;
        long nSteps = 10000;
        // Monson reports data for 0.373773507616 and 0.389566754417
        double density = 0.373773507616;
        String filename = "normal_modes_hexane";
        double fud = 1.0;
        
        boolean graphic = false;

        //parse arguments
        //filename is element 0
        if(args.length > 0){
            filename = args[0];
        }
        if(args.length > 1){
            nSteps = Long.parseLong(args[1]);
        }
        if(args.length > 2){
            density = Double.parseDouble(args[2]);
            if(density == 0.37) {density = 0.373773507616;}
            if(density == 0.40) {density = 0.389566754417;}
        }
        if(args.length > 3){
            fud = Double.parseDouble(args[3]);
        }
        if(args.length > 6){
            xLng = Integer.parseInt(args[3]);
            yLng = Integer.parseInt(args[4]);
            zLng = Integer.parseInt(args[5]);
        }

        
        filename = "nm_hexane";
        nSteps = 100000;
        xLng = 4;
        yLng = 4;
        zLng = 3;
        
        
        
        
        System.out.println("Running " + "3D hexane harmonic simulation, " +
                "measuring hard sphere energy");
        System.out.println(xLng*yLng*zLng + " molecules in a " + xLng + ", " + 
                yLng + ", " + zLng + " arrangement at density " + density );
        System.out.println("Harmonic fudge: " + fud);
        System.out.println(nSteps + " steps");
        
        
        //spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        SimHarmonicHexane sim = new SimHarmonicHexane(Space3D.getInstance(), 
                density, xLng, yLng, zLng, filename, fud);

        //Add hard potentials for FEP calculation.  With de novo sampling
        // potential is not otherwise used.
        Potential p2 = new P2HardSphere(sim.getSpace(), 1.0, true);
        PotentialMaster potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(p2, new AtomType[]
               {((SpeciesHexane)sim.getSpecies(0)).getLeafType(),
                ((SpeciesHexane)sim.getSpecies(0)).getLeafType()});
        /* Yes, Andrew, I know that the above specifies only the first thing in 
         * the array, but this is a hexane simulation, so this should not be an 
         * issue.
         */
        
        if(potentialMaster instanceof PotentialMasterList){
            double L = Math.pow(0.4165783882178116, 1.0/3.0);
            double neighborRange = L / Math.sqrt(2.0);
            ((PotentialMasterList)potentialMaster).setRange(neighborRange);
            //find neighbors now; don't hook up NeighborListManager because
            //the neighbors won't change
            ((PotentialMasterList)potentialMaster).getNeighborManager(sim.box).reset();
        }
        
        //meters for FEP calculations
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, sim.box);
        BoltzmannProcessor bp = new BoltzmannProcessor();
        DataPump pump = new DataPump(meterPE, bp);
        AccumulatorAverage avgBoltzmann = new AccumulatorAverageFixed(1);
        bp.setDataSink(avgBoltzmann);
        avgBoltzmann.setPushInterval(5);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        
        //set up things for determining energy of harmonic system
        //read and set up wave vectors
        
        //GRAPHIC
        if(graphic){
            //meter for harmonic system energy, sent to direct and to 
            // boltzmann average
            MeterHarmonicEnergy harmonicEnergy = new 
                MeterHarmonicEnergy(sim.coordinateDefinition, sim.normalModes);
            DataFork harmonicFork = new DataFork();
            AccumulatorAverage harmonicAvg = new AccumulatorAverageFixed(5);
            DataPump pumpHarmonic = new DataPump(harmonicEnergy, harmonicFork);
            harmonicFork.addDataSink(harmonicAvg);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pumpHarmonic));
            
            //set up measurement of S matrix, to check that configurations are
            // generated as expected
            MeterNormalMode meterNormalMode = new MeterNormalMode();
            meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
            WaveVectorFactory waveVectorFactory = 
                sim.normalModes.getWaveVectorFactory();
            meterNormalMode.setWaveVectorFactory(waveVectorFactory);
            meterNormalMode.setBox(sim.box);
            
            //graphic simulation - set up window
            SimulationGraphic simG = new SimulationGraphic(sim, APP_NAME);
            ArrayList dataStreamPumps = simG.getController().getDataStreamPumps();
            dataStreamPumps.add(pump);
            dataStreamPumps.add(pumpHarmonic);
            
            DisplayTextBoxesCAE boxesPE = new DisplayTextBoxesCAE();
            boxesPE.setAccumulator(avgBoltzmann);
            boxesPE.setPrecision(6);
            simG.add(boxesPE);
            
            DisplayTextBoxesCAE harmonicBoxes = new DisplayTextBoxesCAE();
            harmonicBoxes.setAccumulator(harmonicAvg);
            simG.add(harmonicBoxes);
            
            simG.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
            simG.makeAndDisplayFrame(APP_NAME);
        } else {
        //UNGRAPHIC
            //not graphic, so simulation is run batch & S data is written to file
            sim.activityIntegrate.setMaxSteps(nSteps);
            
            sim.getController().actionPerformed();
            
            DataGroup boltzmannData = (DataGroup)avgBoltzmann.getData();
            double pNotOverlap = ((DataDouble) boltzmannData.getData(AccumulatorAverage.AVERAGE.index)).x;
            double pError = ((DataDouble) boltzmannData.getData(AccumulatorAverage.ERROR.index)).x;
            
            System.out.println("avg HS Boltzmann factor "
                    + pNotOverlap + " +/- " + pError);
            System.out.println("free energy contribution "
                    + (-Math.log(pNotOverlap)) + " +/- " + (pError/pNotOverlap));
            System.out.println("free energy contribution per molecule "
                    + (-Math.log(pNotOverlap)/xLng/yLng/zLng) + " +/- " 
                    + (pError/pNotOverlap)/xLng/yLng/zLng);
        }
            
        
        
        
        PrimitiveHexane primitive = (PrimitiveHexane)sim.lattice.getPrimitive();
         
        
        
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");

        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(nSteps);

        sim.getController().actionPerformed();
                
        System.out.println(filename + " is finished!");
            
    }


}

