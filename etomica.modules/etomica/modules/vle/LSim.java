package etomica.modules.vle;

import etomica.api.IBox;
import etomica.api.ISpecies;

import etomica.action.Action;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterPressureByVolumeChange;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedBox;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;
import etomica.util.HistoryCollapsingAverage;

public class LSim extends Simulation {

    public final IBox boxLiquid;
    public final ISpecies species;
    public final IntegratorMC integratorLiquid;
    public final PotentialMaster potentialMaster;
    public final ActivityIntegrate activityIntegrate;
    public final P2SoftSphericalTruncatedBox potential;
    
    public LSim(int numAtoms, double temperature, double density, Space _space) {
        super(_space);

        double initBoxSize = Math.pow(numAtoms/density, (1.0/3.0));
        
        species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        boxLiquid = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize), space);
        addBox(boxLiquid);
        boxLiquid.setNMolecules(species, numAtoms);
        boxLiquid.setDensity(density);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(), space);
        config.initializeCoordinates(boxLiquid);
        
        potentialMaster = new PotentialMaster(space);
        P2LennardJones p2LJ = new P2LennardJones(space, 1.0, 1.0);
        potential = new P2SoftSphericalTruncatedBox(p2LJ);
        potential.setBox(boxLiquid);
        potentialMaster.addPotential(potential, new ISpecies[]{species, species});
        
        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature);
        integratorLiquid.setBox(boxLiquid);
        MCMoveAtom atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
        atomMove.setStepSize(0.3);
        
        activityIntegrate = new ActivityIntegrate(integratorLiquid);
        getController().addAction(activityIntegrate);
        
        BoxImposePbc pbc = new BoxImposePbc(boxLiquid, space);
        integratorLiquid.addIntervalAction(pbc);
        integratorLiquid.setActionInterval(pbc, 1000);
    }
    
    public static void main(String[] args) {
        Space sp = Space3D.getInstance();
        final LSim sim = new LSim(864, 1.0, 0.70081, sp);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "L", 50, sp);
        simGraphic.getDisplayBox(sim.boxLiquid).setPixelUnit(new Pixel(15));
        MeterPressure meterPressureLiquid = new MeterPressure(sim.getSpace());
        meterPressureLiquid.setIntegrator(sim.integratorLiquid);
        DataFork fork = new DataFork();
        DataPump pumpPressureLiquid = new DataPump(meterPressureLiquid, fork);
        simGraphic.getController().getDataStreamPumps().add(pumpPressureLiquid);
        AccumulatorAverage avgPressureLiquid = new AccumulatorAverageCollapsing();
        avgPressureLiquid.setPushInterval(10);
        fork.addDataSink(avgPressureLiquid);
        AccumulatorHistory historyPressureLiquid = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integratorLiquid);
        historyPressureLiquid.setTimeDataSource(stepCounter);
        fork.addDataSink(historyPressureLiquid);
        historyPressureLiquid.setPushInterval(1);

        DisplayPlot plotHistoryPressure = new DisplayPlot();
        plotHistoryPressure.setLabel("Pressure History");

        DisplayTextBoxesCAE displayPressureLiquid = new DisplayTextBoxesCAE();
        displayPressureLiquid.setLabel("Liquid pressure");
        displayPressureLiquid.setAccumulator(avgPressureLiquid);
        simGraphic.add(displayPressureLiquid);
        historyPressureLiquid.addDataSink(plotHistoryPressure.getDataSet().makeDataSink());
        plotHistoryPressure.setLegend(new DataTag[]{meterPressureLiquid.getTag()}, "Liquid");
        sim.integratorLiquid.addIntervalAction(pumpPressureLiquid);
        sim.integratorLiquid.setActionInterval(pumpPressureLiquid, 500);
        simGraphic.makeAndDisplayFrame();
        simGraphic.add(plotHistoryPressure);
        
        ModifierGeneral modifier = new ModifierGeneral(sim.potential, "truncationFactor");
        DeviceSlider cutoffSlider = new DeviceSlider(sim.getController(), modifier);
        cutoffSlider.setPrecision(2);
        cutoffSlider.setMinimum(0.0);
        cutoffSlider.setMaximum(0.5);
        cutoffSlider.setNMajor(10);
        cutoffSlider.setValue(0.49);
        cutoffSlider.setLabel("cutoff");
        simGraphic.add(cutoffSlider);
    }
    
    public static class NoGUI {
        public static void main(String[] args) {
            int nA = 500;
            double temperature = 1.0;
            double d = 0.70081;
            long steps = 1000000;
            double truncFactor = 0.40;
            if (args.length > 0) {
                nA = Integer.parseInt(args[0]);
                temperature = Double.parseDouble(args[1]);
                d = Double.parseDouble(args[2]);
                steps = Long.parseLong(args[3])*1000;
                truncFactor = Double.parseDouble(args[4]);
            }
            final int numAtoms = nA;
            final double density = d;
            System.out.println("N="+numAtoms);
            System.out.println("T="+temperature);
            System.out.println("rho="+density);
            System.out.println(steps+" steps");
            Space sp = Space3D.getInstance();
            final LSim sim = new LSim(numAtoms, temperature, density, sp);
            sim.potential.setTruncationFactor(truncFactor);
            System.out.println("cutoff="+sim.potential.getRange()+" ("+truncFactor+")");
            sim.activityIntegrate.setMaxSteps(steps);
            sim.getController().actionPerformed();
            System.out.println("equilibration finished");

            sim.getController().reset();

            MeterPressure meterPressureLiquid = new MeterPressure(sim.getSpace());
            meterPressureLiquid.setIntegrator(sim.integratorLiquid);
            final AccumulatorAverage avgPressureLiquid = new AccumulatorAverageCollapsing();
            DataPump pumpPressureLiquid = new DataPump(meterPressureLiquid, avgPressureLiquid);
            sim.integratorLiquid.addIntervalAction(pumpPressureLiquid);
            sim.integratorLiquid.setActionInterval(pumpPressureLiquid, 1000);
            
            final MeterPressureByVolumeChange meterPressureLiquidByDV = new MeterPressureByVolumeChange(sim.getSpace(), new boolean[]{true,true,true});
            meterPressureLiquidByDV.setIntegrator(sim.integratorLiquid);
            DataSplitter splitter = new DataSplitter();
            int n = meterPressureLiquidByDV.getDataInfo().getLength();
            final AccumulatorAverageCollapsing[] avgPressureLiquidByDV = new AccumulatorAverageCollapsing[n]; 
            DataPump pumpPressureLiquidByDV = new DataPump(meterPressureLiquidByDV, splitter);
            sim.integratorLiquid.addIntervalAction(pumpPressureLiquidByDV);
            sim.integratorLiquid.setActionInterval(pumpPressureLiquidByDV, 1000);
            for (int i=0; i<n; i++) {
                avgPressureLiquidByDV[i] = new AccumulatorAverageCollapsing();
                splitter.setDataSink(i, avgPressureLiquidByDV[i]);
            }
            
            Action writeAction = new Action() {
                public void actionPerformed() {
                    double pressure = ((DataGroup)avgPressureLiquid.getData()).getData(AccumulatorAverage.StatType.AVERAGE.index).getValue(0);
                    double error = ((DataGroup)avgPressureLiquid.getData()).getData(AccumulatorAverage.StatType.ERROR.index).getValue(0);
                    System.out.println("P="+pressure+" +/- "+error);
                    
                    double[] vScale = ((DataDoubleArray)meterPressureLiquidByDV.getScalingDataSource().getData()).getData();
                    double volume = numAtoms / density;
                    for (int i=0; i<avgPressureLiquidByDV.length; i++) {
                        double pressureByDV = ((DataGroup)avgPressureLiquidByDV[i].getData()).getData(AccumulatorAverage.StatType.AVERAGE.index).getValue(0);
                        double errorByDV = ((DataGroup)avgPressureLiquidByDV[i].getData()).getData(AccumulatorAverage.StatType.ERROR.index).getValue(0);
                        pressure = Math.log(pressureByDV)/((vScale[i]-1)*volume);
                        error = errorByDV/pressureByDV/(Math.abs(vScale[i]-1)*volume);
                        System.out.println(i+" "+pressure+" "+error);
                    }
                }
            };
            sim.integratorLiquid.addIntervalAction(writeAction);
            sim.integratorLiquid.setActionInterval(writeAction, (int)(steps/10));
            sim.activityIntegrate.setMaxSteps(steps);

            sim.getController().actionPerformed();
            
            double pressure = ((DataGroup)avgPressureLiquid.getData()).getData(AccumulatorAverage.StatType.AVERAGE.index).getValue(0);
            double error = ((DataGroup)avgPressureLiquid.getData()).getData(AccumulatorAverage.StatType.ERROR.index).getValue(0);
            System.out.println("final P="+pressure+" +/- "+error);
        }
    }
}
