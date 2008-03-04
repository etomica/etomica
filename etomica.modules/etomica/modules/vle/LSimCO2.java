package etomica.modules.vle;

import etomica.action.Action;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataTag;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2LJQ;
import etomica.potential.P2SoftTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.HistoryCollapsingAverage;

public class LSimCO2 extends Simulation {

    public final IBox boxLiquid;
    public final SpeciesSpheresRotating species;
    public final IntegratorMC integratorLiquid;
    public final PotentialMaster potentialMaster;
    public final ActivityIntegrate activityIntegrate;
    public final P2SoftTruncated potential;
    public final P2LJQ p2LJ;
    
    public LSimCO2(int numAtoms, double temperature, double density, double cutoffFac) {
        super(Space3D.getInstance());

        int initNumMolecules = 300;
        double sigma = 3.82;
        double epsilon = Kelvin.UNIT.toSim(182.9);
        double moment = Debye.UNIT.toSim(5.1);
        moment *= moment;

        double initBoxSize = Math.pow(initNumMolecules/density, (1.0/3.0));
        
        species = new SpeciesSpheresRotating(this);
        getSpeciesManager().addSpecies(species);
        ((AtomTypeSphere)species.getLeafType()).setDiameter(sigma);

        boxLiquid = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize), space);
        addBox(boxLiquid);
        boxLiquid.setNMolecules(species, numAtoms);
        boxLiquid.setDensity(density);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(), space);
        config.initializeCoordinates(boxLiquid);
        
        potentialMaster = new PotentialMaster(space);
        p2LJ = new P2LJQ(space, sigma, epsilon, moment);
        p2LJ.setTemperature(temperature);
        potential = new P2SoftTruncated(p2LJ, cutoffFac*sigma);
        potential.setBox(boxLiquid);
        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});
        
        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature);
        integratorLiquid.setBox(boxLiquid);
        MCMoveAtom atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
        atomMove.setStepSize(0.3);
        MCMoveRotate rotateMove = new MCMoveRotate(potentialMaster, random);
        integratorLiquid.getMoveManager().addMCMove(rotateMove);
        
        activityIntegrate = new ActivityIntegrate(integratorLiquid);
        getController().addAction(activityIntegrate);
        
        BoxImposePbc pbc = new BoxImposePbc(boxLiquid, space);
        integratorLiquid.addIntervalAction(pbc);
        integratorLiquid.setActionInterval(pbc, 1000);
    }
    
    public static void main(String[] args) {
        final LSimCO2 sim = new LSimCO2(256, Kelvin.UNIT.toSim(273.15), 0.004, 3.0);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "L", 50, sim.space);
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
        
        ModifierGeneral modifier = new ModifierGeneral(sim.potential, "truncationRadius");
        DeviceSlider cutoffSlider = new DeviceSlider(sim.getController(), modifier);
        cutoffSlider.setPrecision(2);
        cutoffSlider.setMinimum(0.0);
        cutoffSlider.setMaximum(sim.boxLiquid.getBoundary().getDimensions().x(1)*0.5);
        cutoffSlider.setNMajor(10);
        cutoffSlider.setValue(3);
        cutoffSlider.setLabel("cutoff");
        simGraphic.add(cutoffSlider);
    }
    
    public static class NoGUI {
        public static void main(String[] args) {
            int nA = 256;
            double temperature = 273.15;
            double d = 0.0006;
            long steps = 10000;
            if (args.length > 0) {
                nA = Integer.parseInt(args[0]);
                temperature = Double.parseDouble(args[1]);
                d = Double.parseDouble(args[2]);
                steps = Long.parseLong(args[3])*1000;
            }
            temperature = Kelvin.UNIT.toSim(temperature);
            final int numAtoms = nA;
            final double density = d;
            System.out.println("N="+numAtoms);
            System.out.println("T="+temperature);
            System.out.println("rho="+density);
            System.out.println(steps+" steps");
            final LSimCO2 sim = new LSimCO2(numAtoms, temperature, density, 3.0);
            System.out.println("boxL = "+sim.boxLiquid.getBoundary().getDimensions().x(0));
            sim.activityIntegrate.setMaxSteps(steps);
            sim.getController().actionPerformed();
            System.out.println("equilibration finished");

            sim.getController().reset();

            final int nCutoff = 5;
            final double cutoffStep = 1;
            final MeterPressure2[][] pressureMeters = new MeterPressure2[nCutoff][2];
            final AccumulatorAverage[][] pressureAvg = new AccumulatorAverage[nCutoff][2];
            for (int i=0; i<nCutoff; i++) {
                PotentialMaster potentialMaster = new PotentialMaster(sim.getSpace());
                P2SoftTruncated potential = new P2SoftTruncated(sim.p2LJ, (i*cutoffStep+2.5)*sim.p2LJ.getSigma());
                potential.setBox(sim.boxLiquid);
                potentialMaster.addPotential(potential, new AtomType[]{sim.species.getLeafType(), sim.species.getLeafType()});
                
                pressureMeters[i][0] = new MeterPressure2(sim.getSpace());
                pressureMeters[i][0].setIncludeLrc(true);
                pressureMeters[i][0].setBox(sim.boxLiquid);
                pressureMeters[i][0].setPotentialMaster(potentialMaster);
                pressureAvg[i][0] = new AccumulatorAverageCollapsing();
                DataPump pump = new DataPump(pressureMeters[i][0], pressureAvg[i][0]);
                sim.integratorLiquid.addIntervalAction(pump);
                sim.integratorLiquid.setActionInterval(pump, 1000);
            
                pressureMeters[i][1] = new MeterPressure2(sim.getSpace());
                pressureMeters[i][1].setIncludeLrc(false);
                pressureMeters[i][1].setBox(sim.boxLiquid);
                pressureMeters[i][1].setPotentialMaster(potentialMaster);
                pressureAvg[i][1] = new AccumulatorAverageCollapsing();
                pump = new DataPump(pressureMeters[i][1], pressureAvg[i][1]);
                sim.integratorLiquid.addIntervalAction(pump);
                sim.integratorLiquid.setActionInterval(pump, 1000);
            }
//            final MeterPressureByVolumeChange meterPressureLiquidByDV = new MeterPressureByVolumeChange(sim.getSpace(), new boolean[]{true,true,true});
//            meterPressureLiquidByDV.setIntegrator(sim.integratorLiquid);
//            DataSplitter splitter = new DataSplitter();
//            int n = meterPressureLiquidByDV.getDataInfo().getLength();
//            final AccumulatorAverageCollapsing[] avgPressureLiquidByDV = new AccumulatorAverageCollapsing[n]; 
//            DataPump pumpPressureLiquidByDV = new DataPump(meterPressureLiquidByDV, splitter);
//            sim.integratorLiquid.addIntervalAction(pumpPressureLiquidByDV);
//            sim.integratorLiquid.setActionInterval(pumpPressureLiquidByDV, 1000);
//            for (int i=0; i<n; i++) {
//                avgPressureLiquidByDV[i] = new AccumulatorAverageCollapsing();
//                splitter.setDataSink(i, avgPressureLiquidByDV[i]);
//            }
            
            Action writeAction = new Action() {
                public void actionPerformed() {
                    for (int i=0; i<nCutoff; i++) {
                        double pressureLRC = ((DataGroup)pressureAvg[i][0].getData()).getData(AccumulatorAverage.StatType.AVERAGE.index).getValue(0);
                        double pressureNoLRC = ((DataGroup)pressureAvg[i][1].getData()).getData(AccumulatorAverage.StatType.AVERAGE.index).getValue(0);
//                        double errorLRC = ((DataGroup)pressureAvg[i][0].getData()).getData(AccumulatorAverage.StatType.ERROR.index).getValue(0);
                        System.out.println((i*cutoffStep+2.5)+" "+pressureLRC+" "+pressureNoLRC);
                    }
                    
//                    double[] vScale = ((DataDoubleArray)meterPressureLiquidByDV.getScalingDataSource().getData()).getData();
//                    double volume = numAtoms / density;
//                    for (int i=0; i<avgPressureLiquidByDV.length; i++) {
//                        double pressureByDV = ((DataGroup)avgPressureLiquidByDV[i].getData()).getData(AccumulatorAverage.StatType.AVERAGE.index).getValue(0);
//                        double errorByDV = ((DataGroup)avgPressureLiquidByDV[i].getData()).getData(AccumulatorAverage.StatType.ERROR.index).getValue(0);
//                        pressure = Math.log(pressureByDV)/((vScale[i]-1)*volume);
//                        error = errorByDV/pressureByDV/(Math.abs(vScale[i]-1)*volume);
//                        System.out.println(i+" "+pressure+" "+error);
//                    }
                }
            };
            sim.integratorLiquid.addIntervalAction(writeAction);
            sim.integratorLiquid.setActionInterval(writeAction, (int)(steps/10));
            sim.activityIntegrate.setMaxSteps(steps);

            sim.getController().actionPerformed();
        }
    }
}
