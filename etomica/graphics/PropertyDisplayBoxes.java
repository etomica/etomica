/*
 * Created on Feb 10, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.graphics;

import java.awt.Component;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.Constants;
import etomica.data.AccumulatorAverage;
import etomica.graphics.DisplayBox.LabelType;
import etomica.units.Unit;

/**
 * @author kenbenjamin
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class PropertyDisplayBoxes extends Display{

	public AccumulatorAverage accumulatorAverage;
	public DisplayBox currentBox, averageBox, errorBox;
	public JPanel panelParentGroup;
	
    protected Constants.CompassDirection labelPosition = Constants.NORTH;
    protected LabelType labelType;

	
	JLabel jLabelPanelParentGroup = new JLabel();
	/**
	 * 
	 */
	public PropertyDisplayBoxes() {
		super();
		currentBox = new DisplayBox();
		averageBox = new DisplayBox();
		errorBox = new DisplayBox();
		jLabelPanelParentGroup = new JLabel();
        panelParentGroup = new JPanel(new java.awt.BorderLayout());
        panelParentGroup.add(currentBox.graphic(), java.awt.BorderLayout.WEST);
        panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.CENTER);
        panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
        setLabel("Label");
        setLabelType(DisplayBox.STRING);
        currentBox.setLabel("Current");
        averageBox.setLabel("Average");
        errorBox.setLabel("Error");

        
        // TODO Auto-generated constructor stub
	}

    /**
     * Specifies the datasourcethat generates the displayed value.
     */
    public void setAccumulator(AccumulatorAverage accumulatorAverage) {
        this.accumulatorAverage = accumulatorAverage;
        currentBox.setDataSource(accumulatorAverage.makeDataSourceMostRecent());
        averageBox.setDataSource(accumulatorAverage.makeDataSourceAverage());
        errorBox.setDataSource(accumulatorAverage.makeDataSourceError());

    }
    
    /**
     * Accessor method for the datsource that generates the displayed value.
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


    public void doUpdate() {
    		currentBox.doUpdate();
    		averageBox.doUpdate();
    		errorBox.doUpdate();
    }

    public void setUnit(Unit unit) {
    		currentBox.setUnit(unit);
    		averageBox.setUnit(unit);
    		errorBox.setUnit(unit);
    }
    
}
