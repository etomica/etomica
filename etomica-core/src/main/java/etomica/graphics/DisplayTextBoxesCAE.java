/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.data.AccumulatorAverage;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.graphics.DisplayTextBox.LabelType;
import etomica.units.Unit;
import etomica.units.dimensions.Null;
import etomica.units.systems.UnitSystem;
import etomica.util.Constants;

import javax.swing.*;
import java.awt.*;

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
    public DisplayTextBox currentBox, averageBox, errorBox, corBox;
	public JPanel panelParentGroup;
    protected boolean doShowCurrent = true, doShowCorrelation = false;
	
    protected Constants.CompassDirection labelPosition = Constants.CompassDirection.NORTH;
    protected LabelType labelType;

	JLabel jLabelPanelParentGroup = new JLabel();
    
    public DisplayTextBoxesCAE() {
		super();
		currentBox = new DisplayTextBox();
        currentBox.setShowUnit(false);
		currentBox.setLabel("Current");
		averageBox = new DisplayTextBox();
        averageBox.setShowUnit(false);
		averageBox.setLabel("Average");
		errorBox = new DisplayTextBox();
        errorBox.setShowUnit(false);
		errorBox.setLabel("Error");
        corBox = new DisplayTextBox();
        corBox.setLabel("Correlation");
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
        if (newDoShowCurrent) setDoShowCorrelation(false);
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

    public void setDoShowCorrelation(boolean newDoShowCorrelation) {
        if (newDoShowCorrelation) setDoShowCurrent(false);
        if (newDoShowCorrelation == doShowCorrelation) return;
        doShowCorrelation = newDoShowCorrelation;
        if (doShowCorrelation) {
            panelParentGroup.remove(averageBox.graphic());
            panelParentGroup.remove(errorBox.graphic());
            panelParentGroup.add(averageBox.graphic(), BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), BorderLayout.CENTER);
            panelParentGroup.add(corBox.graphic(), BorderLayout.EAST);
        } else {
            panelParentGroup.remove(averageBox.graphic());
            panelParentGroup.remove(errorBox.graphic());
            panelParentGroup.remove(corBox.graphic());
            panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.CENTER);
            panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
        }
    }

    public void putDataInfo(IDataInfo dataInfo) {
        if(!(dataInfo instanceof DataInfoGroup)) {
            // we need a data group so we can unpack it and send the pieces on to the actual DisplayTextBoxes
            throw new IllegalArgumentException("DisplayBoxesCAE strangely is being given something other than a DataGroup");
        }
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
                    accumulatorAverage.MOST_RECENT, accumulatorAverage.AVERAGE, accumulatorAverage.ERROR, accumulatorAverage.BLOCK_CORRELATION});
        }
    }

    /**
     * Accessor method for the data source that generates the displayed value.
     */
    public AccumulatorAverage getAccumulator() {
        return accumulatorAverage;
    }

    public Component graphic() {
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
        corBox.setPrecision(3);
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
        DataGroup dg = (DataGroup) data;
        currentBox.putData(dg.getData(0));
        averageBox.putData(dg.getData(1));
        errorBox.putData(dg.getData(2));
        if (dg.getNData() > 3) {
            corBox.putData(((DataGroup) data).getData(3));
        }
    }

    /**
     * Sets the units for the displayed quantities.
     * Does not change the units label for the display.
     */
    public void setUnit(Unit unit) {
        currentBox.setUnit(unit);
        averageBox.setUnit(unit);
        errorBox.setUnit(unit);
    }

    public Unit getUnit() {
        return currentBox.getUnit();
    }
}
