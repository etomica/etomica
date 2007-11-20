package etomica.graphics;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.border.TitledBorder;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

import etomica.action.Action;
import etomica.data.Data;
import etomica.data.DataPump;
import etomica.data.DataSourceFunction;
import etomica.data.DataSourcePoints;
import etomica.modifier.ModifierGeneral;
import etomica.units.Length;
import etomica.units.Null;
import etomica.util.Function;

public class DevicePlotPoints {

	private static final int X_DIM = 0;
	private static final int Y_DIM = 1;
	public static final int MIN_X = 0;
	public static final int MAX_X = 1;
	public static final int MIN_Y = 2;
	public static final int MAX_Y = 3;
	private final int TOTAL_X_SAMPLES = 20;
    private final double[] scaleInit = new double[] {-50.0, 50.0, -250.0, 250.0};

    private double[] scaleMins = new double[] {-50.0, 0.0, -250.0, 0.0};
    private double[] scaleMaxs = new double[] {0.0, 50.0, 0.0, 250.0};

    private JPanel controlPanel;
    private JPanel parameterPanel;
    private JPanel resizePanel;
    private JPanel scalePanel;

    private DeviceTableModelGeneric tableModel;
    private DisplayPlot plot;
    private DeviceTable table;

	private DeviceSlider[] plotSizeSliders = new DeviceSlider[4];
	private JRadioButton buttonAuto;
    private JRadioButton buttonUser;
    private DeviceSlider funcSlider[];
    private Action updateAction;
    private GridBagConstraints vertGBC;

    private int numFunctions;
    private Function[] functions;
    private DataSourceFunction[] dataSourceFuncs;
    private DataPump funcPumps[];
    private ModifierGeneral[] mods;
    private String[] funcParmLabels;

	public DevicePlotPoints(Function[] fncts, String[] funcNames) {
		this(null, fncts, funcNames, true);
	}

	public DevicePlotPoints(String[] sliderLabels,
			Function[] fncts, String[] funcNames, boolean verticalParameters) {

		numFunctions = fncts.length;
		funcParmLabels = sliderLabels;
		int numSliders = 0;
		if(sliderLabels != null) {
		    numSliders = sliderLabels.length;
		}

		// GBC that is reused throughout
		vertGBC = new GridBagConstraints();
    	vertGBC.gridx = 0;
    	vertGBC.gridy = GridBagConstraints.RELATIVE;
        vertGBC.insets = new java.awt.Insets(3,1,3,1);

        // Graphics panels
        controlPanel = new JPanel();
        controlPanel.setLayout(new GridBagLayout());
        resizePanel = new JPanel();
        resizePanel.setLayout(new BorderLayout());

		// function parameter adjustment sliders
	    JPanel sliderPanel[] = new JPanel[numSliders];
	    funcSlider = new DeviceSlider[numSliders];
	    mods = new ModifierGeneral[numSliders];

		for(int slide = 0; slide < numSliders; slide++) {
			sliderPanel[slide] = new JPanel();
			sliderPanel[slide].setBorder(new TitledBorder(null, sliderLabels[slide],
					                     TitledBorder.CENTER, TitledBorder.TOP));
			funcSlider[slide] = new DeviceSlider(null);
			funcSlider[slide].setPrecision(1);
			funcSlider[slide].setMinimum(0);
			funcSlider[slide].setMaximum(10);
			funcSlider[slide].setNMajor(5);
			funcSlider[slide].setValue(5);
			funcSlider[slide].setShowValues(true);
			funcSlider[slide].setEditValues(false);
			mods[slide] = new ModifierGeneral(funcSlider[slide], "value");
			funcSlider[slide].setModifier(mods[slide]);
			sliderPanel[slide].add(funcSlider[slide].graphic());
        }

		//
		// Data point table
		//
		tableModel = new DeviceTableModelGeneric(null, new String[] {"X", "Y"});
		table = new DeviceTable(tableModel);
		table.setPreferredSize(200, 200);
		table.setSize(200, 200);
		table.setTitle("Data Points");
        tableModel.addTableModelListener(new TableChangeListener());

        // Deleting a point
        DeviceButton delButton = new DeviceButton(null); 
        delButton.setLabel("Delete Point(s)");
        ((JPanel)table.graphic()).add(delButton.graphic(), vertGBC);

	    Action deletePointAction = new Action() {
	        public void actionPerformed() {
	        	
	        	int[] selRows = table.getSelectedRows();
	        	while(selRows.length > 0) {
	        		if(selRows[0] == (tableModel.getRowCount()-1)) {
	        			break;
	        		}
	        		else {
	        			tableModel.deleteRow(selRows[0]);
	        			selRows = table.getSelectedRows();
	        		}
	        	}
	            plot.getPlot().repaint();
	        }
	    };
        delButton.setAction(deletePointAction);

        // Add table to the "control panel"
        controlPanel.add(table.graphic(), vertGBC);
        
        // Add parameter adjust scrollbars to a scrolled window which
        // is added to the "control panel"
        JPanel sPanel = new JPanel();
        sPanel.setLayout(new GridBagLayout());
        GridBagConstraints horizGBC = new GridBagConstraints();
        horizGBC.gridx = GridBagConstraints.RELATIVE;
        horizGBC.gridy = 0;
		for(int slide = 0; slide < numSliders; slide++) {
			sPanel.add(sliderPanel[slide], verticalParameters ? vertGBC : horizGBC);
		}
        parameterPanel = new JPanel();
        if (verticalParameters) {
            JScrollPane scrollPane = new JScrollPane(sPanel);
            scrollPane.setPreferredSize(new java.awt.Dimension(250, 320));
            parameterPanel.add(scrollPane);
            controlPanel.add(parameterPanel, vertGBC);
        }
        else {
            parameterPanel.setLayout(new GridBagLayout());
            parameterPanel.add(sPanel);
        }
		parameterPanel.setBorder(new TitledBorder(null, "Function Parameter Adjustment",
                TitledBorder.CENTER, TitledBorder.TOP));

		//
		// Plot axis adjustment sliders
		//
		scalePanel = new JPanel();
		scalePanel.setLayout(new GridBagLayout());
		scalePanel.setBorder(new TitledBorder(null, "Plot Scaling", TitledBorder.CENTER, TitledBorder.TOP));

		JPanel scaleTypePanel = new JPanel();
		scaleTypePanel.setBorder(new TitledBorder(null, "Scaling Type", TitledBorder.CENTER, TitledBorder.TOP));
        ButtonGroup scaleGroup = new ButtonGroup();
        buttonAuto = new JRadioButton("Auto Scale");
        buttonUser = new JRadioButton("User Scale");
        scaleGroup.add(buttonAuto);
        scaleGroup.add(buttonUser);
        buttonUser.setSelected(true);
        scaleTypePanel.add(buttonAuto);
        scaleTypePanel.add(buttonUser);

        JPanel[] scalePanels = new JPanel[4];
        String[] scaleTitles = new String[] {"Minimum", "Maximum", "Minimum", "Maximum"};
        ModifierGeneral[] scaleMods = new ModifierGeneral[4];

        for(int slide = MIN_X; slide <= MAX_Y; slide++) {
	        
			scalePanels[slide] = new JPanel();
			scalePanels[slide].setBorder(new TitledBorder(null, scaleTitles[slide], TitledBorder.CENTER, TitledBorder.TOP));
			plotSizeSliders[slide] = new DeviceSlider(null);
			plotSizeSliders[slide].setPrecision(0);
			plotSizeSliders[slide].setMinimum(scaleMins[slide]);
			plotSizeSliders[slide].setMaximum(scaleMaxs[slide]);
			plotSizeSliders[slide].setNMajor(5);
			plotSizeSliders[slide].setValue(scaleInit[slide]);
			scaleMods[slide] = new ModifierGeneral(plotSizeSliders[slide], "value");
			plotSizeSliders[slide].setModifier(scaleMods[slide]);
            scalePanels[slide].add(plotSizeSliders[slide].graphic());
        }

		JPanel scaleXAxisPanel = new JPanel();
		scaleXAxisPanel.setBorder(new TitledBorder(null, "X Axis Scaling", TitledBorder.CENTER, TitledBorder.TOP));
		scaleXAxisPanel.add(scalePanels[MIN_X]);
		scaleXAxisPanel.add(scalePanels[MAX_X]);

		JPanel scaleYAxisPanel = new JPanel();
		scaleYAxisPanel.setBorder(new TitledBorder(null, "Y Axis Scaling", TitledBorder.CENTER, TitledBorder.TOP));
		scaleYAxisPanel.add(scalePanels[MIN_Y]);
		scaleYAxisPanel.add(scalePanels[MAX_Y]);

		scalePanel.add(scaleTypePanel, vertGBC);
		scalePanel.add(scaleXAxisPanel, vertGBC);
		scalePanel.add(scaleYAxisPanel, vertGBC);

		resizePanel.add(scalePanel, BorderLayout.SOUTH);

		// The Plot
        plot = new DisplayPlot();
        plot.getPlot().setTitle("Function Display");
        plot.setSize(650, 400);

		// Initialize functions displayed on plot
		functions = new Function[numFunctions];
		dataSourceFuncs = new DataSourceFunction[numFunctions];
		funcPumps = new DataPump[numFunctions];
		for(int f = 0; f < numFunctions; f++) {
			functions[f] = fncts[f];
			dataSourceFuncs[f] = new DataSourceFunction(funcNames[f],Null.DIMENSION,functions[f],
					                                    TOTAL_X_SAMPLES,"x",Length.DIMENSION);
			funcPumps[f] = new DataPump(dataSourceFuncs[f], plot.getDataSet().makeDataSink());
		}

		// point display on the plot
		final DataSourcePoints dspts = new DataSourcePoints("Independent Points", Null.DIMENSION);
	    final DataPump ptPump = new DataPump(dspts, plot.getDataSet().makeDataSink());
        plot.getPlot().setMarksStyle("dots", numFunctions);
        plot.getPlot().setConnected(false, numFunctions);

	    //
	    // Update action to pipe any data changes and redraw plot
	    //
	    updateAction = new Action() {
	        public void actionPerformed() {

	            // sniff out min and max y values from the functions.
	            double maxY = Double.MIN_VALUE;
	            double minY = Double.MAX_VALUE;
	            for(int f = 0; f < numFunctions; f++) {
	                dataSourceFuncs[f].getXSource().setXMax(plotSizeSliders[MAX_X].getValue());
	                dataSourceFuncs[f].getXSource().setXMin(plotSizeSliders[MIN_X].getValue());
	                dataSourceFuncs[f].update();
	                Data yFunc = dataSourceFuncs[f].getData();
	                for (int i=0; i<yFunc.getLength(); i++) {
	                    if (yFunc.getValue(i) > maxY) {
	                        maxY = yFunc.getValue(i);
	                    }
                        if (yFunc.getValue(i) < minY) {
                            minY = yFunc.getValue(i);
                        }
	                }

		            funcPumps[f].actionPerformed();
	        	}
	            dspts.update(getPoints(X_DIM), getPoints(Y_DIM));
	            // we could also sniff min and max y here
	            ptPump.actionPerformed();
	            // don't auto-scale in X.  this means that entered data outside
	            // the range won't show up.
                plot.getPlot().setXRange(plotSizeSliders[MIN_X].getValue(),
                        plotSizeSliders[MAX_X].getValue());
                
	            if(buttonUser.isSelected() == true) {
	                plot.getPlot().setYRange(plotSizeSliders[MIN_Y].getValue(),
	            		                     plotSizeSliders[MAX_Y].getValue());
	            }
	            else {
                    plot.getPlot().setYRange(minY, maxY);
	            }
	        }
	    };

        for(int slide = 0; slide < numSliders; slide++) {
            funcSlider[slide].setPostAction(updateAction);
        }

        plotSizeSliders[MIN_X].setPostAction(updateAction);
        plotSizeSliders[MAX_X].setPostAction(updateAction);
        plotSizeSliders[MIN_Y].setPostAction(updateAction);
        plotSizeSliders[MAX_Y].setPostAction(updateAction);

        //
        // Update on scale type selection
        //
        buttonAuto.addActionListener(new ToggleButtonListener());
        buttonUser.addActionListener(new ToggleButtonListener());


        plot.getDataSet().setUpdatingOnAnyChange(true);

        // Update axis min/max values
        updateAction.actionPerformed();

	}

    /**
     * Returns the top level panel that the function plot components sit on.
     */
    public DisplayPlot getDisplayPlot() {
        return plot;
    }

	/**
	 * Returns the top level panel that the function plot components sit on.
	 */
	public Component graphic() {
		return plot.graphic();
	}

	/**
	 * Returns the top level panel that the points table and (if parameters are
	 * vertically-oriented) function parameter slider components sit on.
	 */
	public Component controlGraphic() {
		return controlPanel;
	}

    /**
     * Returns the panel that contains the parameter sliders.  You probably
     * don't want this unless you have horizontal layout for the parameters.
     */
    public Component parameterGraphic() {
        return parameterPanel;
    }

	/**
	 * Returns the top level panel that the plot scale resizing components sit on.
	 */
	public Component resizeGraphic() {
		return resizePanel;
	}

	/**
	 * Instructs the plot to redisplay any plotted functions with the most
	 * recent data.  Any points contained in the data point input table will be
	 * displayed (if the data point input table is displayed).  The plot axis
	 * limits are updated. 
	 */
	public void refresh() {
		updateAction.actionPerformed();
	}

	/**
	 * Returns the plot scale slider for the given parameter (values can be
	 * DevicePlotPoints.MIN_X, MAX_X, MIN_Y or MAX_Y.
	 */
	public DeviceSlider getPlotSizeSlider(int minMaxXY) {
	    return plotSizeSliders[minMaxXY];
	}
	
	/**
	 * Sets the minimum X axis limit allowed on the function plot.  The slider that
	 * sets the minimum X axis limit will have a range of this value to zero.
	 * @param min  The minimum X axis value
	 */
	public void setMinimumXScale(double min) {
        if (scaleMaxs[MAX_X] < min) {
            setMaximumXScale(min);
        }
		scaleMins[MIN_X] = min;
		plotSizeSliders[MIN_X].setMinimum(scaleMins[MIN_X]);
	}

    /**
     * Sets the minimum value for the maximum X slider.
     */
    public void setMinimumMaxXScale(double minMax) {
        if (scaleMins[MIN_X] > minMax) {
            setMinimumXScale(minMax);
        }
        if (scaleMaxs[MAX_X] < minMax) {
            setMaximumXScale(minMax);
        }
        scaleMins[MAX_X] = minMax;
        plotSizeSliders[MAX_X].setMinimum(scaleMins[MAX_X]);
    }

    /**
     * Sets the maxmimum value for the minimum X slider.
     */
    public void setMaximumMinXScale(double maxMin) {
        if (scaleMaxs[MAX_X] < maxMin) {
            setMaximumXScale(maxMin);
        }
        if (scaleMins[MIN_X] > maxMin) {
            setMinimumXScale(maxMin);
        }
        scaleMaxs[MIN_X] = maxMin;
        plotSizeSliders[MIN_X].setMaximum(scaleMaxs[MIN_X]);
    }

	/**
	 * Sets the maximum X axis limit allowed on the function plot.  The slider that
	 * sets the maximum X axis limit will have a range of zero to this value.
	 * @param max  The maximum X axis value
	 */
	public void setMaximumXScale(double max) {
	    if (scaleMins[MIN_X] > max) {
	        setMinimumXScale(max);
	    }
		scaleMaxs[MAX_X] = max;
		plotSizeSliders[MAX_X].setMaximum(scaleMaxs[MAX_X]);
	}

	/**
	 * Sets the minimum Y axis limit allowed on the function plot.  The slider that
	 * sets the minimum Y axis limit will have a range of this value to zero.
	 * @param min  The minimum Y axis value
	 */
	public void setMinimumYScale(double min) {
		scaleMins[MIN_Y] = min;
		plotSizeSliders[MIN_Y].setMinimum(scaleMins[MIN_Y]);
	}

	/**
	 * Sets the maximum Y axis limit allowed on the function plot.  The slider that
	 * sets the maximum Y axis limit will have a range of zero to this value.
	 * @param max  The maximum Y axis value
	 */
	public void setMaximumYScale(double max) {
		scaleMaxs[MAX_Y] = max;
		plotSizeSliders[MAX_Y].setMaximum(scaleMaxs[MAX_Y]);
	}

	/**
	 * Sets the visibility of the plot axis limit sliders.  Default is true(show).
	 * @param show show/unshow boolean
	 */
	public void showScale(boolean show) {
		if(show == false) {
		    resizePanel.remove(scalePanel);
		}
		else {
			resizePanel.add(scalePanel, BorderLayout.SOUTH);
		}
	}

	/**
	 * Sets the visibility of the data point input table.  Default is true(show).
	 * @param show show/unshow boolean
	 */
	public void showPointInput(boolean show) {
		if(show == false) {
	        controlPanel.remove(table.graphic());
		}
		else {
			controlPanel.add(table.graphic(), vertGBC);
		}
	}

	/**
	 * Returns the value of the function parameter slider.
	 * @param desc label of parameter setting slider.  The label was passed
	 * into the ctor.
	 * @return slider setting of the parameter value.
	 */
	public double getParameterValue(String desc) {
		double parmValue = 0.0;
		for(int i = 0; i < funcParmLabels.length; i++) {
		    if(desc.compareTo(funcParmLabels[i]) == 0) {
		    		parmValue = mods[i].getValue();
		    		break;
		    }
		}
		return parmValue;
	}

	/**
	 * Sets the minimum and maximum allowable values for a function parameter.
	 * @param desc label of parameter setting slider.  The label was passed
	 * into the ctor.
	 * @param min Minimum allowable value for the function parameter.
	 * @param max Maxminimum allowable value for the function parameter.
	 */
	public void setParameterLimits(String desc, double min, double max) {

		for(int i = 0; i < funcParmLabels.length; i++) {
		    if(desc.compareTo(funcParmLabels[i]) == 0) {
		    	funcSlider[i].setMinimum(min);
		    	funcSlider[i].setMaximum(max);
		    	if(funcSlider[i].getValue() < min) {
		    		funcSlider[i].setValue(min);
		    	}
		    	else if(funcSlider[i].getValue() > max) {
		    		funcSlider[i].setValue(max);
		    	}
		    	break;
		    }
		}
	}

	private double[] getPoints(int column) {

        double[] points = null;
        int nonBlankRowCount = 0;
        int numRows = tableModel.getRowCount();

		for(int row = 0; row < numRows; row++) {
			if(((String)tableModel.getValueAt(row, column)).compareTo("") != 0) {
                nonBlankRowCount++;
			}
		}

		points = new double[nonBlankRowCount];
		for(int row = 0; row < numRows; row++) {
			if(((String)tableModel.getValueAt(row, column)).compareTo("") != 0) {
				points[row] = Double.valueOf(((String)tableModel.getValueAt(row, column))).doubleValue();
			}
		}

		return points;
	}
	
	/**
	 * Returns the slider for the given parameter.
	 */
	public DeviceSlider getSlider(String desc) {
        for(int i = 0; i < funcParmLabels.length; i++) {
            if(desc.compareTo(funcParmLabels[i]) == 0) {
                return funcSlider[i];
            }
        }
        return null;
	}
	
	public void setAutoScale(boolean isAutoScale) {
	    buttonAuto.setSelected(isAutoScale);
	    new ToggleButtonListener().actionPerformed(null);
	    updateAction.actionPerformed();
	}

    private class ToggleButtonListener implements ActionListener {
    	public void actionPerformed(ActionEvent e) {
            if(buttonAuto.isSelected()) {
            	plotSizeSliders[MIN_X].getSlider().setEnabled(false);
            	plotSizeSliders[MAX_X].getSlider().setEnabled(false);
            	plotSizeSliders[MIN_Y].getSlider().setEnabled(false);
            	plotSizeSliders[MAX_Y].getSlider().setEnabled(false);
            }
            else {
            	plotSizeSliders[MIN_X].getSlider().setEnabled(true);
            	plotSizeSliders[MAX_X].getSlider().setEnabled(true);
            	plotSizeSliders[MIN_Y].getSlider().setEnabled(true);
            	plotSizeSliders[MAX_Y].getSlider().setEnabled(true);
            }		
            updateAction.actionPerformed();
        }
    }

    // TODO : Hmmm, no order gaurenteed between this listener and other listeners...
    // Like the one in DeviceTableModelGeneric
    private class TableChangeListener implements TableModelListener {
    	public void tableChanged(TableModelEvent e) {
    		// If a row changed and is a complete entry(both x and y values)
    		// plot the points on the display.
    		if(e.getType() == TableModelEvent.UPDATE) {
        		Object blank = "";

        		if(tableModel.getValueAt(e.getFirstRow(), 0).equals(blank) == false &&
        		   tableModel.getValueAt(e.getFirstRow(), 1).equals(blank) == false) {
                    updateAction.actionPerformed();
                }
    		}
    		// If a row is removed, remove it from the display.
        	else if(e.getType() == TableModelEvent.DELETE) {
                updateAction.actionPerformed();
        	}
    	}
    } // end class TableChangeListener


}
