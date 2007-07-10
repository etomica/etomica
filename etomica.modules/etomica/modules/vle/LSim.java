package etomica.modules.vle;

import etomica.action.Action;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataTag;
import etomica.data.meter.MeterPressure;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.IPotential;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.HistoryCollapsingAverage;

public class LSim extends Simulation {

    public final Box boxLiquid;
    public final Species species;
    public final IntegratorMC integratorLiquid;
    public final PotentialMasterCell potentialMaster;
    public final ActivityIntegrate activityIntegrate;
    public final P2SoftSphericalTruncated potential;
    protected double sigma;
    protected double temperature;
    protected double epsilon;
    protected double density;
    
    public LSim() {
        super(Space3D.getInstance());
        int initNumMolecules = 2000;
        sigma = 1.0;
        temperature = 1.0;
        epsilon = 1.0;
        density = 0.3;
        double cutoff = 6*sigma;

        double initBoxSize = Math.pow(initNumMolecules/density, (1.0/3.0));
        
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);

        boxLiquid = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize));
        addBox(boxLiquid);
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxLiquid.setDensity(0.70081);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc());
        config.initializeCoordinates(boxLiquid);
        
        potentialMaster = new PotentialMasterCell(this, cutoff);
        P2LennardJones p2LJ = new P2LennardJones(space, sigma, epsilon);
        potential = new P2SoftSphericalTruncated(p2LJ, cutoff);
        potentialMaster.addPotential(potential, new Species[]{species, species});
        potentialMaster.setCellRange(3);
        
        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature);
        integratorLiquid.setBox(boxLiquid);
        MCMoveAtom atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        
        activityIntegrate = new ActivityIntegrate(integratorLiquid, false, false);
        getController().addAction(activityIntegrate);

        integratorLiquid.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(boxLiquid).makeMCMoveListener());
        potentialMaster.getNbrCellManager(boxLiquid).assignCellAll();
    }
    
    public static void main(String[] args) {
        final LSim sim = new LSim();
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "L", 50);
        MeterPressure meterPressureLiquid = new MeterPressure(sim.getSpace());
        meterPressureLiquid.setIntegrator(sim.integratorLiquid);
        DataFork fork = new DataFork();
        DataPump pumpPressureLiquid = new DataPump(meterPressureLiquid, fork);
        simGraphic.getController().getDataStreamPumps().add(pumpPressureLiquid);
        AccumulatorAverage avgPressureLiquid = new AccumulatorAverage(40);
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
        cutoffSlider.setMinimum(0);
        cutoffSlider.setMaximum(6);
        cutoffSlider.setLabel("cutoff");
        cutoffSlider.setPostAction(new Action(){
            public void actionPerformed() {
                sim.potentialMaster.setRange(sim.potential.getTruncationRadius());
                sim.potentialMaster.reset();
                sim.potentialMaster.getNbrCellManager(sim.boxLiquid).assignCellAll();
            }
        });
        simGraphic.add(cutoffSlider);
    }
}
