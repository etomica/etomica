/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.Angstrom;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;

import javax.swing.*;
import java.util.ArrayList;

/**
 * Graphic UI for Droplet module.  Design by Ludwig Nitsche.
 *
 * @author Andrew Schultz
 */
public class DropletAtomicGraphicFasterer extends SimulationGraphic {

    private final static String APP_NAME = "Droplet";
    private final static int REPAINT_INTERVAL = 1;
    protected DropletAtomicFasterer sim;

    public DropletAtomicGraphicFasterer(final DropletAtomicFasterer simulation, Space _space) {

        super(simulation, TABBED_PANE, APP_NAME, _space.D() == 2 ? 10 * REPAINT_INTERVAL : REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        this.sim = simulation;

        getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                sim.makeDropShape();
            }
        });

        sim.potentialMaster.getNeighborManager().setDoDownNeighbors(true);
        final AtomTestLiquidAtomicFasterer liquidFilter = new AtomTestLiquidAtomicFasterer(sim.potentialMaster.getNeighborManager(), sim.box);
        sim.potentialMaster.init();
        final ColorSchemeLiquidVapor colorScheme = new ColorSchemeLiquidVapor(liquidFilter);
        colorScheme.setDoResetFilter(true);
        getDisplayBox(sim.box).setColorScheme(colorScheme);
        ((DiameterHashByType) getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), sim.sigma);
        final DeviceButton cutawayButton = new DeviceButton(sim.getController());
        cutawayButton.setLabel("Liquid Atoms");
        cutawayButton.setAction(new IAction() {
            public void actionPerformed() {
                if (filterIsActive) {
                    getDisplayBox(sim.box).setAtomTestDoDisplay(null);
                } else {
                    getDisplayBox(sim.box).setAtomTestDoDisplay(liquidFilter);
                }
                filterIsActive = !filterIsActive;
                colorScheme.setDoResetFilter(!filterIsActive);

                cutawayButton.setLabel(filterIsActive ? "All Atoms" : "Liquid Atoms");
            }

            boolean filterIsActive = false;
        });
        cutawayButton.setPostAction(getPaintAction(sim.box));
        ((JPanel) getController().graphic()).add(cutawayButton.graphic());

        DisplayTimerFasterer displayTimer = new DisplayTimerFasterer(sim.integrator);
        getPanel().controlPanel.add(displayTimer.graphic(), SimulationPanel.getVertGBC());

        //add meter and display for current kinetic temperature

        DeviceSlider radiusSlider = new DeviceSlider(sim.getController());
        radiusSlider.setPrecision(1);
        radiusSlider.setMinimum(0.2);
        radiusSlider.setMaximum(0.8);
        radiusSlider.setNMajor(4);
        radiusSlider.setShowBorder(true);
        radiusSlider.setLabel("Droplet Size");
        radiusSlider.setModifier(new Modifier() {

            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            public String getLabel() {
                return "Droplet Radius (" + Angstrom.UNIT.symbol() + ")";
            }

            public double getValue() {
                return sim.dropRadius;
            }

            public void setValue(double newValue) {
                sim.dropRadius = newValue;
                sim.makeDropShape();

                if (sim.integrator.getStepCount() > 0) {
                    sim.potentialMaster.init();
                    sim.integrator.reset();
                }

                getDisplayBox(sim.box).repaint();
            }
        });
        //  add(radiusSlider);

        DeviceSlider squeezeSlider = new DeviceSlider(sim.getController());
        squeezeSlider.setShowBorder(true);
        squeezeSlider.setModifier(new ModifierGeneral(sim.p1Smash, "g"));
        squeezeSlider.setPrecision(1);
        squeezeSlider.setMinimum(0);
        squeezeSlider.setMaximum(10);
        squeezeSlider.setNMajor(4);
        squeezeSlider.setShowValues(true);
        squeezeSlider.setLabel("Squeezing force (" + Angstrom.UNIT.symbol() + "/ps^2)");
        add(squeezeSlider);
        float dropDiameter = (float) (0.5 * sim.dropRadius * sim.getBox(0).getBoundary().getBoxSize().getX(0));
        final EllipseDisplayAction ellipseDisplayAction = new EllipseDisplayAction(this, dropDiameter);
        squeezeSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                ellipseDisplayAction.displayEllipse(sim.p1Smash.g);
            }
        });

//        nSlider = new DeviceNSelector(sim.getController());
//        nSlider.setSpecies(sim.species);
//        nSlider.setBox(sim.box);
//        nSlider.setMinimum(0);
//        nSlider.setMaximum(2000);
//        nSlider.setLabel("Number of Atoms");
//        nSlider.setShowBorder(true);
//        nSlider.setShowValues(true);
//        // add a listener to adjust the thermostat interval for different
//        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
//        // don't need as much thermostating.
//        final ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
//        nSlider.setPostAction(new IAction() {
//            public void actionPerformed() {
//                config.initializeCoordinates(sim.box);
//                ((PotentialMasterList)sim.integrator.getPotentialMaster()).getNeighborManager(sim.box).reset();
//                sim.integrator.reset();
//
//                getController().getSimRestart().actionPerformed();
//                getDisplayBox(sim.box).repaint();
//            }
//            
//        });

//        JPanel systemPanel = new JPanel(new GridBagLayout());
//        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        //************* Lay out components ****************//

//        JTabbedPane tabbedPane = new JTabbedPane();
//        tabbedPane.add("System", systemPanel);
//        getPanel().controlPanel.add(tabbedPane, vertGBC);
//        JPanel dropletPanel = new JPanel(new GridBagLayout());
//        numMoleculesPanel.add(nSlider.graphic(), vertGBC);
//        dropletPanel.add(radiusSlider.graphic());
//        tabbedPane.add("Droplet", dropletPanel);
//        JPanel potentialPanel = new JPanel(new GridBagLayout());
//        tabbedPane.add("Surfactant potential", potentialPanel);
//        getPanel().controlPanel.add(tabbedPane);

        DataSourceCountTimeFasterer timeCounter = new DataSourceCountTimeFasterer(sim.integrator);

        MeterPotentialEnergyFromIntegratorFasterer meterPE = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
        DataProcessor foo = new DataProcessor() {

            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                return inputDataInfo;
            }

            protected IData processData(IData inputData) {
                ((DataDouble) inputData).x /= sim.box.getMoleculeList().size() * Kelvin.UNIT.toSim(118);
                return inputData;
            }
        };
        DataFork peFork = new DataFork();
        DataPump pePump = new DataPump(meterPE, foo);
        foo.setDataSink(peFork);
        dataStreamPumps.add(pePump);
        IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
        sim.integrator.getEventManager().addListener(pePumpListener);
        pePumpListener.setInterval(10);
        AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        peFork.addDataSink(peHistory);
        peHistory.setTimeDataSource(timeCounter);
//        AccumulatorAverageCollapsing peAvg = new AccumulatorAverageCollapsing();
//        peFork.addDataSink(peAvg);
//        DisplayTextBoxesCAE peBox = new DisplayTextBoxesCAE();
//        peBox.setAccumulator(peAvg);
//        peAvg.setPushInterval(1);
//        add(peBox);

//        MeterKineticEnergy meterKE = new MeterKineticEnergy();
//        meterKE.setBox(sim.box);
//        DataFork keFork = new DataFork();
//        DataPump kePump = new DataPump(meterKE, keFork);
//        dataStreamPumps.add(kePump);
//        sim.integrator.addIntervalAction(kePump);
//        sim.integrator.setActionInterval(kePump, 10);
//        AccumulatorHistory keHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//        keFork.addDataSink(keHistory);
//        keHistory.setTimeDataSource(timeCounter);

//        MeterEnergy meterE = new MeterEnergy(sim.potentialMaster, sim.box);
//        DataFork eFork = new DataFork();
//        DataPump ePump = new DataPump(meterE, eFork);
//        dataStreamPumps.add(ePump);
//        sim.integrator.addIntervalAction(ePump);
//        sim.integrator.setActionInterval(ePump, 10);
//        AccumulatorHistory eHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//        eFork.addDataSink(eHistory);
//        eHistory.setTimeDataSource(timeCounter);

        DisplayPlotXChart ePlot = new DisplayPlotXChart();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
//        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
//        eHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.getPlot().setYLabel("Potential Energy");
        ePlot.setDoLegend(false);
        ePlot.setLabel("Energy");
        add(ePlot);

        MeterDeformation meterDeformation = new MeterDeformation(space);
        meterDeformation.setBox(sim.box);
        meterDeformation.setFilter(liquidFilter);
        DataSplitter dSplitter = new DataSplitter();
        DataPump dPump = new DataPump(meterDeformation, dSplitter);
        dataStreamPumps.add(dPump);
        IntegratorListenerAction dPumpListener = new IntegratorListenerAction(dPump);
        sim.integrator.getEventManager().addListener(dPumpListener);
        dPumpListener.setInterval(10);

        DataFork dFork = new DataFork();
        dSplitter.setDataSink(1, dFork);
        AccumulatorHistory dHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        dFork.addDataSink(dHistory);
        dHistory.setTimeDataSource(timeCounter);

        DisplayPlotXChart deformationPlot = new DisplayPlotXChart();
        dHistory.setDataSink(deformationPlot.getDataSet().makeDataSink());
        deformationPlot.setLabel("Deformation");
        deformationPlot.getPlot().setYLabel("Deformation");
        deformationPlot.setDoLegend(false);
        add(deformationPlot);

        DataFork rFork = new DataFork();
        dSplitter.setDataSink(0, rFork);
        AccumulatorHistory rHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        rFork.addDataSink(rHistory);
        rHistory.setTimeDataSource(timeCounter);

        DisplayPlotXChart radiusPlot = new DisplayPlotXChart();
        rHistory.setDataSink(radiusPlot.getDataSet().makeDataSink());
        radiusPlot.setLabel("Radius");
        radiusPlot.getPlot().setYLabel("Radius (" + Angstrom.UNIT.symbol() + ")");
        radiusPlot.setDoLegend(false);
        add(radiusPlot);

//        
//        DataSourceScalar liquidDensity = new DataSourceScalar("Liquid Density", new DimensionRatio(Quantity.DIMENSION, Volume.DIMENSION)) {
//            private static final long serialVersionUID = 1L;
//
//            public double getDataAsScalar() {
//                CellLattice lattice = sim.potentialMaster.getNbrCellManager(sim.box).getLattice();
//                Cell cell = (Cell)lattice.site(v);
//                
//                double[] size = lattice.getCellSize();
//                return cell.occupants().getAtomCount() / (size[0]*size[1]*size[2]);
//            }
//            
//            protected final Vector v = space.makeVector();
//        };
//        
//        DataFork ldFork = new DataFork();
//        DataPump ldPump = new DataPump(liquidDensity, ldFork);
//        dataStreamPumps.add(ldPump);
//        sim.integrator.addIntervalAction(ldPump);
//        sim.integrator.setActionInterval(ldPump, 10);
//        AccumulatorHistory ldHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//        ldFork.addDataSink(ldHistory);
//        ldHistory.setTimeDataSource(timeCounter);
//        
//        DataSourceScalar vaporDensity = new DataSourceVaporDensity("Vapor Density", new DimensionRatio(Quantity.DIMENSION, Volume.DIMENSION));
//        
//        DataFork vdFork = new DataFork();
//        DataPump vdPump = new DataPump(vaporDensity, vdFork);
//        dataStreamPumps.add(vdPump);
//        sim.integrator.addIntervalAction(vdPump);
//        sim.integrator.setActionInterval(vdPump, 10);
//        AccumulatorHistory vdHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//        vdFork.addDataSink(vdHistory);
//        vdHistory.setTimeDataSource(timeCounter);
//        
//        DisplayPlotXChart dPlot = new DisplayPlotXChart();
//        ldHistory.setDataSink(dPlot.getDataSet().makeDataSink());
//        vdHistory.setDataSink(dPlot.getDataSet().makeDataSink());
//        dPlot.setLabel("Density");
//        add(dPlot);
    }

    public static void main(String[] args) {
        Space sp = null;
        if (args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    sp = Space3D.getInstance();
                } else {
                    sp = Space2D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        } else {
            sp = Space3D.getInstance();
        }

        DropletAtomicFasterer sim = new DropletAtomicFasterer();
        DropletAtomicGraphicFasterer swGraphic = new DropletAtomicGraphicFasterer(sim, sp);
        swGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(2));
        SimulationGraphic.makeAndDisplayFrame
                (swGraphic.getPanel(), APP_NAME);
    }
}


