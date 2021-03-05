/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


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
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.units.Pixel;

import java.util.ArrayList;

/**
 * Simulation to sample harmonic potential
 */
public class SimHarmonic extends Simulation {

	private static final String APP_NAME = "Sim Harmonic";
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;

    public Box box;
    public Boundary boundary;
    public SpeciesGeneral species;
    public NormalModes normalModes;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    public SimHarmonic(Space _space, int numAtoms, double density, String filename, double harmonicFudge) {
        super(_space);

        int D = space.D();

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0 / density);
            boundary = new BoundaryRectangularPeriodic(space, numAtoms / density);
            nCells = new int[]{numAtoms};
        } else {
            primitive = new PrimitiveFcc(space, 1);
            double v = primitive.unitCell().getVolume();
            primitive.scaleSize(Math.pow(v * density, -1.0 / 3.0));
            int n = (int) Math.round(Math.pow(numAtoms, 1.0 / 3.0));
            nCells = new int[]{n, n, n};
            boundary = new BoundaryDeformableLattice(primitive, nCells);
        }
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(this, null, box);

        this.getController().addActivity(new ActivityIntegrate(integrator));

        MCMoveHarmonic move = new MCMoveHarmonic(getRandom());
        integrator.getMoveManager().addMCMove(move);


        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, space);
        coordinateDefinition.initializeCoordinates(nCells);

        if (D == 1) {
            normalModes = new NormalModes1DHR(boundary, numAtoms);
        } else {
            normalModes = new NormalModesFromFile(filename, D);
        }
        normalModes.setTemperature(1.0);
        normalModes.setHarmonicFudge(harmonicFudge);

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        move.setOmegaSquared(normalModes.getOmegaSquared());
        move.setEigenVectors(normalModes.getEigenvectors());
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinition);
        move.setTemperature(1.0);

        move.setBox(box);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        //set up simulation parameters
        int D = 1;
        int nA = 27;
        double density = 1.3;
        double harmonicFudge = 1.0;
        long steps = 800000;
        if (D == 1) {
            nA = 3;
            density = 0.5;
        }
        boolean graphic = false;
        String filename = "normal_modes3D_27_130";
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            steps = Long.parseLong(args[2]);
        }
        if (args.length > 3) {
            nA = Integer.parseInt(args[3]);
        }
        if (args.length > 4) {
            harmonicFudge = Double.parseDouble(args[4]);
        }
        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" harmonic simulation, measuring hard sphere energy");
        System.out.println(nA+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println(steps+" MC steps");

        //construct simulation
        SimHarmonic sim = new SimHarmonic(Space.getInstance(D), nA, density, filename, harmonicFudge);

        //add hard potentials for FEP calculations.  With de novo sampling potential is not otherwise used.
        Potential2 p2 = new P2HardSphere(sim.getSpace(), 1.0, true);
        if (D == 1) {
            p2 = new P2XOrder(sim.getSpace(), (Potential2HardSpherical)p2);
        }
        PotentialMaster potentialMaster = (D == 1 ? new PotentialMasterList(sim, sim.space) : new PotentialMasterMonatomic(sim.getSpeciesManager()));
        potentialMaster.addPotential(p2, new AtomType[]{sim.species.getLeafType(), sim.species.getLeafType()});

        if (potentialMaster instanceof PotentialMasterList) {
            double neighborRange;
            if (D == 1) {
                neighborRange = 1.01 / density;
            }
            else {
                //FCC
                double L = Math.pow(0.26*density, 1.0/3.0);
                neighborRange = L / Math.sqrt(2.0);
            }
            ((PotentialMasterList)potentialMaster).setRange(neighborRange);
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMaster).getNeighborManager(sim.box).reset();
        }

        //meters for FEP calculations
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, sim.box);
        BoltzmannProcessor bp = new BoltzmannProcessor();
        bp.setTemperature(1);
        DataPump pump = new DataPump(meterPE,bp);
        AccumulatorAverage avgBoltzmann = new AccumulatorAverageFixed(1);
        bp.setDataSink(avgBoltzmann);
        avgBoltzmann.setPushInterval(5);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));

//         MeterMomentumCOM meterCOM = new MeterMomentumCOM(sim.space);
//         MeterPositionCOM meterCOM = new MeterPositionCOM(sim.space);
//         DataSinkConsole console = new DataSinkConsole();
//         DataProcessorFunction filter = new DataProcessorFunction(new Function.Chop());
//         DataPump comPump = new DataPump(meterCOM,filter);
//         filter.setDataSink(console);
//         IntervalActionAdapter comAdapter = new IntervalActionAdapter(comPump);
//         sim.integrator.addListener(comAdapter);
//         meterCOM.setBox(sim.box);

        //set up things for determining energy of harmonic system
        //read and set up wave vectors

        if(graphic){
            //meter for harmonic system energy, sent to direct and to boltzmann average
            MeterHarmonicEnergy harmonicEnergy = new MeterHarmonicEnergy(sim.coordinateDefinition, sim.normalModes);
            DataFork harmonicFork = new DataFork();
            AccumulatorAverage harmonicAvg = new AccumulatorAverageFixed(5);
            DataPump pumpHarmonic = new DataPump(harmonicEnergy, harmonicFork);
            harmonicFork.addDataSink(harmonicAvg);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pumpHarmonic));

            //histogram energy of individual modes
//            MeterHarmonicSingleEnergy harmonicSingleEnergy = new MeterHarmonicSingleEnergy(coordinateDefinitionLeaf, sim.normalModes);
//            harmonicSingleEnergy.setTemperature(1.0);
//            harmonicSingleEnergy.setBox(sim.box);
//    //        DataProcessorFunction harmonicLog = new DataProcessorFunction(new Function.Log());
//            AccumulatorAverage harmonicSingleAvg = new AccumulatorAverage(5);
//            DataHistogram harmonicSingleHistogram = new DataHistogram(new HistogramSimple.Factory(50, new DoubleRange(0, 1)));
//            pump = new DataPump(harmonicSingleEnergy, harmonicSingleHistogram);
//    //        harmonicLog.setDataSink(harmonicSingleHistogram);
//            harmonicSingleHistogram.setDataSink(harmonicSingleAvg);
//            iaa= new IntervalActionAdapter(pump);
//            iaa.setActionInterval(1);
//            sim.integrator.addListener(iaa);

            //set up measurement of S matrix, to check that configurations are generated as expected
            MeterNormalMode meterNormalMode = new MeterNormalMode();
            meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
            WaveVectorFactory waveVectorFactory = sim.normalModes.getWaveVectorFactory();
            meterNormalMode.setWaveVectorFactory(waveVectorFactory);
            meterNormalMode.setBox(sim.box);


            //graphic simulation -- set up window
//            sim.getDefaults().pixelUnit = new Pixel(0.05);
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

//            DisplayPlot harmonicPlot = new DisplayPlot();
//            harmonicPlot.setDoLegend(false);
//            harmonicSingleAvg.addDataSink(harmonicPlot.getDataSet().makeDataSink(), new StatType[]{StatType.AVERAGE});
//            simG.add(harmonicPlot);

            simG.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
            simG.makeAndDisplayFrame(APP_NAME);
        } else {
            //not graphic, so run simulation batch
            //S data is written to file
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

            DataGroup boltzmannData = (DataGroup)avgBoltzmann.getData();
            double pNotOverlap = ((DataDouble) boltzmannData.getData(avgBoltzmann.AVERAGE.index)).x;
            double pError = ((DataDouble) boltzmannData.getData(avgBoltzmann.ERROR.index)).x;

            System.out.println("avg HS Boltzmann factor "+pNotOverlap+" +/- "+pError);

            System.out.println("free energy contribution "+(-Math.log(pNotOverlap))+" +/- "+(pError/pNotOverlap));
            System.out.println("free energy contribution per molecule "+(-Math.log(pNotOverlap)/nA)+" +/- "+(pError/pNotOverlap)/nA);
        }

    }
}
