/*
 * (c) Copyright 1996, KL GROUP INC.  
 * ALL RIGHTS RESERVED 
 *
 * This file is provided for demonstration and educational uses only.
 * Permission to use, copy, modify and distribute this file for
 * any purpose and without fee is hereby granted, provided that the
 * above copyright notice and this permission notice appear in all
 * copies, and that the name of KL Group not be used in advertising
 * or publicity pertaining to this material without the specific,
 * prior written permission of an authorized representative of
 * KL Group.
 *
 * KL GROUP MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES, EX-
 * PRESS OR IMPLIED, WITH RESPECT TO THE SOFTWARE, INCLUDING, BUT
 * NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR ANY PARTICULAR PURPOSE, AND THE WARRANTY AGAINST IN-
 * FRINGEMENT OF PATENTS OR OTHER INTELLECTUAL PROPERTY RIGHTS.  THE
 * SOFTWARE IS PROVIDED "AS IS", AND IN NO EVENT SHALL KL GROUP OR
 * ANY OF ITS AFFILIATES BE LIABLE FOR ANY DAMAGES, INCLUDING ANY
 * LOST PROFITS OR OTHER INCIDENTAL OR CONSEQUENTIAL DAMAGES 
 * RELATING TO THE USE OF THIS SOFTWARE.
 */

//   RCSID -- $RCSfile$ $Revision$
//            $Date$  $Locker$  KL Group Inc.

package simulate;

import java.awt.Event;
import java.util.Observable;
import java.util.Vector;
import jclass.chart.ChartDataModelUpdate;
import jclass.chart.EditableChartable;
import jclass.chart.EventTrigger;
import jclass.chart.JCChart;
import jclass.chart.JCChartUtil;
import jclass.table3.TableDataEvent;
import jclass.table3.TableDataListener;
import jclass.util.JCVector;

/**
 * This view models the energy data of an EnergyModel instance. 
 * The TableDataListener lets us know when one of the models 
 * has been changed by calling dataChanged();
 * We can update the models by setting Energy values
 * in setDataItem()
 */

public class EnergyChartView extends Observable
               implements EditableChartable, TableDataListener {

double y[] = {500, 1000, 1500, 2000, 2500}; 
double y1[] = {0, 0, 0, 0, 0};
double y2[] = {0, 0, 0, 0, 0};
double y3[] = {0, 0, 0, 0, 0};
double y4[] = {0, 0, 0, 0, 0};

EnergyDataModel modelEnergy;

public EnergyChartView() {
    this(null);
}

public EnergyChartView(EnergyDataModel one) {
    if (one == null) {
	    JCVector pressure = new JCVector();
	    JCVector configs = new JCVector();
	    JCVector current = new JCVector();
	    JCVector average = new JCVector();
	    JCVector error = new JCVector();
    	for (int i=0; i<y.length; i++) {
	        configs.addElement(new Double((double)y[i]));
	        configs.addElement(new Double((double)y1[i]));
	        current.addElement(new Double((double)y2[i]));
	        average.addElement(new Double((double)y3[i]));
	        error.addElement(new Double((double)y4[i]));
	    }
        modelEnergy = new EnergyDataModel(pressure, configs, current, average, error);
    }
    else modelEnergy = one;
    modelEnergy.addTableDataListener(this);
}

/**
 * When the model data changes, dataChanged() will
 * be called. Just let our observers(the bar chart)
 * know so that they can re-extract the new values.
 */
public void dataChanged(TableDataEvent ev) {
    //System.out.println("ContribChart.dataChanged "+ev.getSource());
    Object source = ev.getSource();
    if (source instanceof EnergyDataModel) {
	int row = (ev.getRow() > 0) ? ev.getRow() : 0;
	int col = (ev.getColumn() > 0) ? ev.getColumn() : 0;
	setChanged();
	//RESET tells chart to re-read all the data of the data source. 
	notifyObservers(new ChartDataModelUpdate(ChartDataModelUpdate.RESET,
						 row, col));
    }
}

/**
 * Gets the value of the Name property.  The Name property defines
 * the name of the data object
 */
    public java.lang.String getName() {
	return "";
    }

/**
 * Sets a data item in the data source. Since DataSource is
 * an Observable object, setDataItem() will cause a model-view
 * update that will pass the new value back to the Observer.
 * This will allow DataSource-side constraints to be placed on 
 * data. 
 * @param row row index
 * @param column column index
 * @param item new data item, not necessarily the final value 
 * assigned to item (row, column)
 * @return a boolean signifying whether or not the change
 * was accepted
 */
    public boolean setDataItem(int row,int column,java.lang.Object item) {
	    if (!(item instanceof Number)) return false;
	    if (row == EnergyDataModel.PRESSURE) return modelEnergy.setTableDataItem(item, EnergyDataModel.PRESSURE, column);
	    if (row == EnergyDataModel.CONFIGS) return modelEnergy.setTableDataItem(item, EnergyDataModel.CONFIGS, column);
        if (row == EnergyDataModel.CURRENT) return modelEnergy.setTableDataItem(item, EnergyDataModel.CURRENT, column); 
	    if (row == EnergyDataModel.AVERAGE) return modelEnergy.setTableDataItem(item, EnergyDataModel.AVERAGE, column);
        if (row == EnergyDataModel.ERROR) return modelEnergy.setTableDataItem(item, EnergyDataModel.ERROR, column); 
	    return false;
    }

/**
 * Gets the value of the SeriesLabel property. The SeriesLabel
 * property defines the label that is associated with a data
 * series. The SeriesLabel can be unparsed JCString text,
 * whereas SeriesName is a name used for retrieving data
 * series.
 * @param row index of the row for which the label is returned
 * @return row/series label
 */
    public String getSeriesLabel(int row) {
	return this.getSeriesName(row);
    }

/**
 * Gets the value of the SeriesName property. The SeriesName
 * property is a string associated with a row in a data source.
 * It is used by ChartDataView as the name for a data series.
 * @param row index of the row for which the name is returned
 * @return row/series name
 */
    public String getSeriesName(int row) {
	switch (row) {
	case 0:
	    return "Current";
	case 1:
	    return "Average";
	case 2:
	    return "Error";
	default:
	}
	return null;
    }

/**
 * Gets the value of the PointLabels property. The PointLabels
 * property is an indexed property that contains a list of point
 * labels provided by the data source. 
 * @return array of point labels
 */
    public java.lang.String[] getPointLabels() {
	JCVector configs = modelEnergy.getConfigs();
	int size = configs.size();
	String[] labels = new String[size];
	Double conf;
	String label;
	for (int i = 0; i < size; i++) {
	    conf = (Double) configs.elementAt(i);
	    label = JCChartUtil.format(conf.doubleValue(), 0);
	    labels[i] = new String(label);
	}
	return labels;
    }

/**
 * Gets the NumRows property of the data source, which
 * specifies the number of rows of data available.
 * @return number of available rows
 */
    public int getNumRows() {
	return 5;
    }

/**
 * Retrieves an entire row of data items.
 * @param row row index
 * @return array of object-derived data items
 */
    public java.util.Vector getRow(int row) {
	switch (row) {
	case 0:
	    return (Vector) modelEnergy.getPressure();
	case 1:
	    return (Vector) modelEnergy.getConfigs();
	case 2:
	    return (Vector) modelEnergy.getCurrent();
	case 3:
	    return (Vector) modelEnergy.getAverage();
	case 4:
	    return (Vector) modelEnergy.getError();
	default:
	}
	return null;
    }

/**
 * Retrieves a data item from the data source.
 * @param row row index
 * @param column column index
 * @return object-derived data item
 */
    public java.lang.Object getDataItem(int row,int column) {
	Vector line = (Vector) this.getRow(row);
	if (line == null) return null;
	return line.elementAt(column);
    }

/**
 * Is this Chartable object holding data in ARRAY or GENERAL format?
 * @return either ARRAY or GENERAL
 */
    public int getDataInterpretation() {
	return ARRAY;
    }
}
