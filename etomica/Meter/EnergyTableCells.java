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

import jclass.table3.EditableVectorDataSource;
import jclass.table3.JCTblEnum;
import jclass.table3.Table;
import jclass.table3.TableDataEvent;
import jclass.table3.TableDataListener;
import jclass.util.JCVector;

/*
 * A Table cell data structure that provides a view on EnergyDataModel
 * When the model gets updated, our TableDataListener dataChanged()
 * will get called and all views will be notified.
 * The Cells also allow the contribution to be set by passing
 * on the update to the model.
 */

public class EnergyTableCells extends EditableVectorDataSource
              implements TableDataListener {

    EnergyDataModel modelEnergy;

    int PRESSURE = 0;
    int CONFIGS = 1;
    int CURRENT = 2;
    int AVERAGE = 3;
    int ERROR = 4;

    double y[] = {500};


    public EnergyTableCells() {
        this(null);
    }

    public EnergyTableCells(EnergyDataModel one) {
        if (one == null) {
	        JCVector pressure = new JCVector();
	        JCVector configs = new JCVector();
	        JCVector current = new JCVector();
	        JCVector average = new JCVector();
	        JCVector error = new JCVector();
    	    for (int i=0; i<y.length; i++) {
	            configs.addElement(new Double((double)y[i]));
	            current.addElement(new Double(i));
	            average.addElement(new Double(i));
	            error.addElement(new Double(i));
	        }
        modelEnergy = new EnergyDataModel(pressure, configs, current, average, error);
        }
       else {
            modelEnergy = one;
       }
        modelEnergy.addTableDataListener(this);

        setColumnLabel(PRESSURE, "Pressure (bar)");
        setColumnLabel(CONFIGS, "Temperature (K)");
        setColumnLabel(CURRENT, "Density (mol/l)");
        setColumnLabel(AVERAGE,   "Avg. Temp.");
        setColumnLabel(ERROR,  "Avg. Dens.");
        this.rows = modelEnergy.getNumColumns();
        this.columns = modelEnergy.getNumRows();
    }

    //called by model data class EnergyDataModel when its data has changed.
    public void dataChanged(TableDataEvent ev) {
        Object source = ev.getSource();
        if (source instanceof EnergyDataModel) {
	        EnergyDataModel model = (EnergyDataModel)source;
	        JCVector pressure = model.getPressure();
	        JCVector configs = model.getConfigs();
	        JCVector current = model.getCurrent();
	        JCVector average = model.getAverage();
	        JCVector error = model.getError();
            //property "cell" declared in VectorDataSource, setCell is method of
            // EditableVectorDataSource an is called below to assign the value in "cell"
            // in this class as opposed to "cell" in the Data Model.
	        for (int row=0; row < modelEnergy.getNumColumns(); row++) {
		        setCell(row, PRESSURE, pressure.elementAt(row));
		        setCell(row, CONFIGS, configs.elementAt(row));
		        setCell(row, CURRENT, current.elementAt(row));
		        setCell(row, AVERAGE, average.elementAt(row));
		        setCell(row, ERROR, error.elementAt(row));
	        }
	    }
	    /**
	    * from VectorDataSource (TableDataSupport), notifies all TableDataListeners (Data Models)
	    * that the data has changed bye calling their dataChanged methods.
	    **/
        setDataChanged(JCTblEnum.ALL, JCTblEnum.ALL, 0, 0, TableDataEvent.RESET);
        this.rows = modelEnergy.getNumColumns();
    }
 
    /**
    * Sets the cell value in the appropriate Data Model (if more than one)
    */
    public boolean setTableDataItem(Object value, int row, int column) {
        //System.out.println("editedItem("+row+","+column+")="+value);
        boolean result = false;
        if (column == PRESSURE) {
	        result = modelEnergy.setTableDataItem(value, modelEnergy.PRESSURE, row);
        }
        else if (column == CONFIGS) {
	        result = modelEnergy.setTableDataItem(value, modelEnergy.CONFIGS, row);
        }
        else if (column == CURRENT) {
    	    result = modelEnergy.setTableDataItem(value, modelEnergy.CURRENT, row);
        }
        else if (column == AVERAGE) {
    	    result = modelEnergy.setTableDataItem(value, modelEnergy.AVERAGE, row);
        }
        else if (column == ERROR) {
    	    result = modelEnergy.setTableDataItem(value, modelEnergy.ERROR, row);
        }
        return result;
    }
}