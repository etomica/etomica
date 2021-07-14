/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.materialfracture;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterPressureTensorFasterer;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure2D;

import javax.swing.border.TitledBorder;
import java.awt.*;

/**
 * Graphical components for Material Fracture module
 */
public class MaterialFractureGraphicFasterer extends SimulationGraphic {

    public MaterialFractureGraphicFasterer(final MaterialFractureFasterer sim) {
        super(sim, SimulationGraphic.TABBED_PANE, "Material Fracture", 10);

        getDisplayBox(sim.box).setPixelUnit(new Pixel(6));
        ((DiameterHashByType) getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), 3.0);

        final StrainColorScheme strainColor = new StrainColorScheme();
        getDisplayBox(sim.box).setColorScheme(strainColor);
        // 37 is the index of the first atom (on the left) to be colored red
        strainColor.setNumber(37);

        getController().getSimRestart().setConfiguration(sim.config);

        DeviceThermoSlider thermoSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        thermoSlider.setIsothermalButtonsVisibility(false);
        thermoSlider.setMaximum(600);
        add(thermoSlider);

        final MeterStrain meterStrain = new MeterStrain();
        meterStrain.setBox(sim.box);
        meterStrain.setAtomNumber(37);
        DataFork strainFork = new DataFork();
        DataPump strainPump = new DataPump(meterStrain, strainFork);
        getController().getDataStreamPumps().add(strainPump);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(strainPump));
        AccumulatorAverageCollapsing strainAverage = new AccumulatorAverageCollapsing();
        strainFork.addDataSink(strainAverage);

        final MeterStressFasterer meterStress = new MeterStressFasterer(sim);
        meterStress.setBox(sim.box);
        DataFork stressFork = new DataFork();
        DataPump stressPump = new DataPump(meterStress, stressFork);
        getController().getDataStreamPumps().add(stressPump);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(stressPump));
        AccumulatorHistory stressHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(5000));
        stressHistory.setTimeDataSource(meterStrain);
        stressFork.addDataSink(stressHistory);
        AccumulatorAverageCollapsing stressAverage = new AccumulatorAverageCollapsing();
        stressAverage.setPushInterval(10);
        stressFork.addDataSink(stressAverage);

        final DisplayPlotXChart stressHistoryPlot = new DisplayPlotXChart();
        stressHistory.setDataSink(stressHistoryPlot.getDataSet().makeDataSink());
        stressHistoryPlot.setLabel("Stress");
        stressHistoryPlot.setDoDrawLines(new DataTag[]{stressHistory.getTag()}, false);
        stressHistoryPlot.getChart().getStyler().setMarkerSize(4);

        add(stressHistoryPlot);

        DisplayTextBoxesCAE stressDisplay = new DisplayTextBoxesCAE();
        stressDisplay.setAccumulator(stressAverage);
        add(stressDisplay);

        DisplayTextBoxesCAE strainDisplay = new DisplayTextBoxesCAE();
        strainDisplay.setAccumulator(strainAverage);
        add(strainDisplay);

        final MeterPressureTensorFasterer meterPressure = new MeterPressureTensorFasterer(sim.integrator.getPotentialCompute(), sim.box);
        DataProcessor pressureToStress = new DataProcessor() {
            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                return dataInfo;
            }

            protected IData processData(IData inputData) {
                // xx component is the first one
                data.x = -inputData.getValue(0);
                return data;
            }

            protected final IDataInfo dataInfo = new DataDouble.DataInfoDouble("Stress", Pressure2D.DIMENSION);
            protected final DataDouble data = new DataDouble();
        };

        DataPump internalStressPump = new DataPump(meterPressure, pressureToStress);
        getController().getDataStreamPumps().add(internalStressPump);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(internalStressPump));

        DataFork internalStressFork = new DataFork();
        pressureToStress.setDataSink(internalStressFork);
        AccumulatorHistory internalStressHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(5000));
        internalStressHistory.setTimeDataSource(meterStrain);
        internalStressFork.addDataSink(internalStressHistory);
        AccumulatorAverageCollapsing internalStressAverage = new AccumulatorAverageCollapsing();
        internalStressAverage.setPushInterval(10);
        internalStressFork.addDataSink(internalStressAverage);

        internalStressHistory.setDataSink(stressHistoryPlot.getDataSet().makeDataSink());
        stressHistoryPlot.setDoDrawLines(new DataTag[]{internalStressHistory.getTag()}, false);
        stressHistoryPlot.setLegend(new DataTag[]{internalStressHistory.getTag()}, "Internal Stress");

        DisplayTextBoxesCAE internalStressDisplay = new DisplayTextBoxesCAE();
        internalStressDisplay.setLabel("Internal Stress");
        internalStressDisplay.setAccumulator(internalStressAverage);
        add(internalStressDisplay);

        final javax.swing.JTabbedPane integratorPanel = new javax.swing.JTabbedPane();
        integratorPanel.setBorder(new TitledBorder(
                null, "Integrator", TitledBorder.CENTER, TitledBorder.TOP));

        ModifierGeneral stepModifier = new ModifierGeneral(sim.integrator, "timeStep");
        DeviceSlider stepSlider = new DeviceSlider(sim.getController(), stepModifier);
        stepSlider.setUnit(new PrefixedUnit(Prefix.FEMTO, Second.UNIT));
        stepSlider.setMaximum(10);
        stepSlider.setNMajor(5);
        stepSlider.setPrecision(1);
        stepSlider.setShowValues(true);
        integratorPanel.add("Timestep (fs)", stepSlider.graphic());

        Modifier nuModifier = new Modifier() {
            public void setValue(double newValue) {
                sim.integrator.setThermostatInterval((int) (1.0 / newValue));
            }

            public double getValue() {
                return 1.0 / sim.integrator.getThermostatInterval();
            }

            public String getLabel() {
                return "";
            }

            public Dimension getDimension() {
                return Null.DIMENSION;
            }
        };
        DeviceSlider nuSlider = new DeviceSlider(sim.getController(), nuModifier);
        nuSlider.setMaximum(1);
        nuSlider.setPrecision(2);
        nuSlider.setNMajor(5);
        nuSlider.setValue(0.1);
        nuSlider.setShowValues(true);
        integratorPanel.add("Andersen nu", nuSlider.graphic());

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        getPanel().controlPanel.add(integratorPanel, vertGBC);

        final javax.swing.JTabbedPane potentialPanel = new javax.swing.JTabbedPane();
        potentialPanel.setBorder(new TitledBorder(
                null, "Potential", TitledBorder.CENTER, TitledBorder.TOP));

        ModifierGeneral epsilonModifier = new ModifierGeneral(sim.p2LJ, "epsilon");
        DeviceSlider epsilonSlider = new DeviceSlider(sim.getController(), epsilonModifier);
        epsilonSlider.setMaximum(50);
        epsilonSlider.setUnit(new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT));
        epsilonSlider.setNMajor(5);
        epsilonSlider.setShowValues(true);
        potentialPanel.add("Epsilon (kJ/mol)", epsilonSlider.graphic());

        ModifierGeneral springConstantModifier = new ModifierGeneral(sim.p1Tension, "springConstant");
        DeviceSlider springConstantSlider = new DeviceSlider(sim.getController(), springConstantModifier);
        springConstantSlider.setUnit(new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT));
        springConstantSlider.setPrecision(2);
        springConstantSlider.setMaximum(0.3);
        springConstantSlider.setNMajor(3);
        springConstantSlider.setShowValues(true);
        potentialPanel.add("Spring Constant (J/(mol A^2))", springConstantSlider.graphic());

        ModifierGeneral cutoffModifier = new ModifierGeneral(sim.pt, "truncationRadius");
        DeviceSlider cutoffSlider = new DeviceSlider(sim.getController(), cutoffModifier);
        cutoffSlider.setMinimum(2);
        cutoffSlider.setMaximum(7);
        cutoffSlider.setNMajor(5);
        cutoffSlider.setPrecision(1);
        cutoffSlider.setShowValues(true);
        potentialPanel.add("Cutoff (A)", cutoffSlider.graphic());

        getPanel().controlPanel.add(potentialPanel, vertGBC);

        getController().getResetAveragesButton().setPostAction(new IAction() {
            public void actionPerformed() {
                stressHistoryPlot.getPlot().clearData();
            }
        });
    }

    /**
     * /* main method
     */
    public static void main(String[] args) {
        MaterialFractureFasterer sim = new MaterialFractureFasterer();
        MaterialFractureGraphicFasterer simGraphic = new MaterialFractureGraphicFasterer(sim);

        simGraphic.makeAndDisplayFrame();
    }
}
