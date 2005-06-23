package etomica.graphics;

import java.awt.Component;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.Constants;
import etomica.Data;
import etomica.DataSink;
import etomica.data.AccumulatorAverage;
import etomica.data.DataGroup;
import etomica.data.DataPump;
import etomica.data.meter.MeterPressureHard;
import etomica.graphics.DisplayBox.LabelType;
import etomica.integrator.IntervalActionAdapter;
import etomica.simulations.HSMD2D;
import etomica.units.Dimension;
import etomica.units.Unit;

/**
 * Display that presents three boxes with the current value, average, 
 * and error, respectively, given by an AccumulatorAverage.
 * 
 * @author kenbenjamin
 */

/*
 * Created on Feb 10, 2005
 */
public class DisplayBoxesCAE extends Display implements DataSink {

	public AccumulatorAverage accumulatorAverage;
	public DisplayBox currentBox, averageBox, errorBox;
	public JPanel panelParentGroup;
	
    protected Constants.CompassDirection labelPosition = Constants.NORTH;
    protected LabelType labelType;

	JLabel jLabelPanelParentGroup = new JLabel();
    
	public DisplayBoxesCAE() {
		super();
		currentBox = new DisplayBox();
		averageBox = new DisplayBox();
		errorBox = new DisplayBox();
		jLabelPanelParentGroup = new JLabel();
        panelParentGroup = new JPanel(new java.awt.BorderLayout());
        panelParentGroup.add(currentBox.graphic(), java.awt.BorderLayout.WEST);
        panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.CENTER);
        panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
        setLabel("");
        setLabelType(DisplayBox.STRING);
        currentBox.setLabel("Current");
        averageBox.setLabel("Average");
        errorBox.setLabel("Error");
	}

    /**
     * Specifies the accumulator that generates the displayed values.
     * Sets up this display to receive the current, average, and error
     * values from the accumulator.
     */
    public void setAccumulator(AccumulatorAverage accumulatorAverage) {
        this.accumulatorAverage = accumulatorAverage;
        accumulatorAverage.addDataSink(this,new AccumulatorAverage.Type[] {
                AccumulatorAverage.MOST_RECENT,
                AccumulatorAverage.AVERAGE,
                AccumulatorAverage.ERROR});
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
	
    public void setLabel(String s) {
        jLabelPanelParentGroup.setText(s);
        if(labelType == DisplayBox.BORDER) {
            panelParentGroup.setBorder(new javax.swing.border.TitledBorder(s));
        }
        if(labelType == DisplayBox.STRING) setLabelPosition(labelPosition);
    }
    
    public String getLabel() {
        return jLabelPanelParentGroup.getText();
    }

    /**
     * Sets label to the given value if it was not previously set.
     * If setLabel was previously called, this method has no effect.
     * This method is usually invoked automatically when this data
     * sink is attached to a data pipe.
     */
    public void setDefaultLabel(String defaultLabel) {
        if(getLabel() == "") setLabel(defaultLabel);
    }


    public void setLabelPosition(Constants.CompassDirection position) {
        labelPosition = position;
        if(labelType != DisplayBox.STRING) return;
        panelParentGroup.remove(jLabelPanelParentGroup);
        panelParentGroup.add(jLabelPanelParentGroup,position.toString());//toString() returns the corresponding BorderLayout constant
//        support.firePropertyChange("label",oldLabel,label);
        panelParentGroup.revalidate();
        panelParentGroup.repaint();
    }
    
    public void setLabelType(LabelType labelType) {
        this.labelType = labelType;
        if(labelType != DisplayBox.BORDER) panelParentGroup.setBorder(new javax.swing.border.EmptyBorder(2,2,2,2));
        if(labelType != DisplayBox.STRING) panelParentGroup.remove(jLabelPanelParentGroup);
        setLabel(jLabelPanelParentGroup.getText());
    }

    public void putData(Data data) {
        currentBox.putData(((DataGroup)data).getData(0));
        averageBox.putData(((DataGroup)data).getData(1));
        errorBox.putData(((DataGroup)data).getData(2));
    }
    
    public void setDimension(Dimension dimension) {
        if(unit == Unit.UNDEFINED) setUnit(dimension.defaultIOUnit());
    }

    public void setUnit(Unit unit) {
    		currentBox.setUnit(unit);
    		averageBox.setUnit(unit);
    		errorBox.setUnit(unit);
    }
    
    private Unit unit = Unit.UNDEFINED;
    
    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.integrator.setIsothermal(true);
        MeterPressureHard pressureMeter = new MeterPressureHard(sim.space,sim.integrator);
        pressureMeter.setPhase(sim.phase);
        AccumulatorAverage accumulator = new AccumulatorAverage();
        DataPump dataPump = new DataPump(pressureMeter, accumulator);
        new IntervalActionAdapter(dataPump, sim.integrator);
        DisplayBoxesCAE boxes = new DisplayBoxesCAE();
        boxes.setAccumulator(accumulator);
        graphic.add(boxes);
        graphic.makeAndDisplayFrame();
    }

}
