/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Component;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.data.*;
import etomica.data.meter.MeterPressureHard;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.graphics.DisplayTextBox.LabelType;
import etomica.listener.IntegratorListenerAction;
import etomica.simulation.prototypes.HSMD2D;
import etomica.units.dimensions.Null;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;
import etomica.util.Constants;

/**
 * Display that presents three boxes with the current value, average, 
 * and error, respectively, given by an AccumulatorAverage.
 * 
 * @author kenbenjamin
 */

/*
 * Created on Feb 10, 2005
 */
public class DisplayTextBoxesCAE extends Display implements IDataSink {

	public AccumulatorAverage accumulatorAverage;
	public DisplayTextBox currentBox, averageBox, errorBox;
	public JPanel panelParentGroup;
	protected boolean doShowCurrent = true;
	
    protected Constants.CompassDirection labelPosition = Constants.CompassDirection.NORTH;
    protected LabelType labelType;

	JLabel jLabelPanelParentGroup = new JLabel();
    
    public DisplayTextBoxesCAE() {
		super();
		currentBox = new DisplayTextBox();
		currentBox.setLabel("Current");
		averageBox = new DisplayTextBox();
		averageBox.setLabel("Average");
		errorBox = new DisplayTextBox();
		errorBox.setLabel("Error");
		jLabelPanelParentGroup = new JLabel();
        panelParentGroup = new JPanel(new java.awt.BorderLayout());
        panelParentGroup.add(currentBox.graphic(), java.awt.BorderLayout.WEST);
        panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.CENTER);
        panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
        setLabel("");
        setLabelType(LabelType.STRING);
        setPrecision(4);
	}

    public void setDoShowCurrent(boolean newDoShowCurrent) {
        if (newDoShowCurrent == doShowCurrent) return;
        doShowCurrent = newDoShowCurrent;
        if (doShowCurrent) {
            panelParentGroup.add(currentBox.graphic(), java.awt.BorderLayout.WEST);
        }
        else {
            panelParentGroup.remove(currentBox.graphic());
        }
    }

    public boolean getDoShowCurrent() {
        return doShowCurrent;
    }

    public DataPipe getDataCaster(IDataInfo dataInfo) {
        if(!(dataInfo instanceof DataInfoGroup)) {
            throw new IllegalArgumentException("DisplayBoxesCAE strangely is being given something other than a DataGroup");
        }
        return null;
    }

    public void putDataInfo(IDataInfo dataInfo) {
        if(getLabel().equals("")) {
            setLabel(dataInfo.getLabel());
        }
        if(getUnit() == Null.UNIT) {
            setUnit(dataInfo.getDimension().getUnit(UnitSystem.SIM));
        }
    }

    /**
     * Specifies the accumulator that generates the displayed values.
     * Sets up this display to receive the current, average, and error
     * values from the accumulator.
     */
    public void setAccumulator(AccumulatorAverage accumulatorAverage) {
        if (this.accumulatorAverage != null) {
            this.accumulatorAverage.removeDataSink(this);
        }
        this.accumulatorAverage = accumulatorAverage;
        if (accumulatorAverage != null) {
            accumulatorAverage.addDataSink(this,new AccumulatorAverage.StatType[] {
                    accumulatorAverage.MOST_RECENT, accumulatorAverage.AVERAGE, accumulatorAverage.ERROR});
        }
    }

    /**
     * Accessor method for the data source that generates the displayed value.
     */
    public AccumulatorAverage getAccumulator() {
        return accumulatorAverage;
    }

	public Component graphic(Object obj){
        return panelParentGroup;
	}

    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public int getPrecision() {
        return currentBox.getPrecision();
    }

    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public void setPrecision(int n) {
        currentBox.setPrecision(n);
        averageBox.setPrecision(n);
        errorBox.setPrecision(n > 2 ? 2 : n);
    }

    public void setLabel(String s) {
        jLabelPanelParentGroup.setText(s);
        if(labelType == LabelType.BORDER) {
            panelParentGroup.setBorder(BorderFactory.createTitledBorder(s));
            
        }
        if(labelType == LabelType.STRING) setLabelPosition(labelPosition);
    }

    public String getLabel() {
        return jLabelPanelParentGroup.getText();
    }

    public void setLabelPosition(Constants.CompassDirection position) {
        labelPosition = position;
        if(labelType != LabelType.STRING) return;
        panelParentGroup.remove(jLabelPanelParentGroup);
        panelParentGroup.add(jLabelPanelParentGroup,position.toString());//toString() returns the corresponding BorderLayout constant
//        support.firePropertyChange("label",oldLabel,label);
        panelParentGroup.revalidate();
        panelParentGroup.repaint();
    }

    public void setLabelType(LabelType labelType) {
        this.labelType = labelType;
        if(labelType != LabelType.BORDER) panelParentGroup.setBorder(BorderFactory.createEmptyBorder(2,2,2,2));
        if(labelType != LabelType.STRING) panelParentGroup.remove(jLabelPanelParentGroup);
        setLabel(jLabelPanelParentGroup.getText());
    }

    public void putData(IData data) {
        currentBox.putData(((DataGroup)data).getData(0));
        averageBox.putData(((DataGroup)data).getData(1));
        errorBox.putData(((DataGroup)data).getData(2));
    }

    public void setUnit(Unit unit) {
        currentBox.setUnit(unit);
        averageBox.setUnit(unit);
        errorBox.setUnit(unit);
    }

    public Unit getUnit() {
        return currentBox.getUnit();
    }

    public static void main(String[] args) {
    	final String APP_NAME = "Display Boxes CAE";

    	etomica.space.Space sp = etomica.space2d.Space2D.getInstance();
        final HSMD2D sim = new HSMD2D();
        final SimulationGraphic graphic = new SimulationGraphic(sim, APP_NAME, sp, sim.getController());
        sim.integrator.setIsothermal(true);
        MeterPressureHard pressureMeter = new MeterPressureHard(sim.getSpace());
        pressureMeter.setIntegrator(sim.integrator);
        AccumulatorAverageCollapsing accumulator = new AccumulatorAverageCollapsing();
        DataPump dataPump = new DataPump(pressureMeter, accumulator);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(dataPump));
        graphic.getController().getDataStreamPumps().add(dataPump);
        DisplayTextBoxesCAE boxes = new DisplayTextBoxesCAE();
        boxes.setAccumulator(accumulator);
        graphic.add(boxes);
        graphic.getController().getReinitButton().setPostAction(graphic.getPaintAction(sim.box));

        graphic.makeAndDisplayFrame(APP_NAME);
    }
}
