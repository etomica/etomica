/*
 * (c) Copyright 1997, 1998 KL GROUP INC.
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
import jclass.table3.TableDataEvent;
import jclass.util.JCVector;
import java.util.Vector;


/*
 * This class represents the base data in this application:
 * An array of ages
 * An array of contributions
 * An array representing the accumulated value of the contribs
 *
 * When the rate or contributions change, the value array
 * is recalculated and all views on the data are notified.
 *
 * This class is based on EditableVectorDataSource
 * but it could just as easily be based on EditableChartable
 */

public class EnergyDataModel extends EditableVectorDataSource {

    double temp = 300;

    public static int CONFIGS = 1;
    public static int CURRENT = 2;
    public static int AVERAGE = 3;
    public static int ERROR = 4;
    public static int PRESSURE = 0;

    public EnergyDataModel(JCVector pressure, JCVector configs, JCVector current, JCVector average, JCVector error) {
        if (configs == null) {
    	    // Default data
	        double y[] = {500};
	        JCVector sample = new JCVector();
	        JCVector sample1 = new JCVector();
	        JCVector sample2 = new JCVector();
	        JCVector sample3 = new JCVector();
	        JCVector sample4 = new JCVector();
	        for (int i=0; i<y.length; i++) {
	            sample.addUnique(new Double(i));
	            sample1.addUnique(new Double(i));
	            sample2.addUnique(new Double(i));
	            sample3.addUnique(new Double(i));
	            sample4.addUnique(new Double(i));
	        }
	        configs = sample1;
	        current = sample2;
	        average = sample3;
	        error = sample4;
	        pressure = sample;
        }
        columns = configs.size();
        columns = 1;
        //cells property declared in VectorDataSource as a JCVector
        cells.addUnique(new JCVector(pressure));
        cells.addUnique(new JCVector(configs));
        cells.addUnique(new JCVector(current));
        cells.addUnique(new JCVector(average));
        cells.addUnique(new JCVector(error));
        update();
    }

    public EnergyDataModel() {
        this(null, null, null, null, null);
    }

    public JCVector getConfigs() {
        return new JCVector((JCVector)cells.elementAt(CONFIGS));
    }

    public JCVector getCurrent() {
        return new JCVector((JCVector)cells.elementAt(CURRENT));
    }
    public JCVector getAverage() {
        return new JCVector((JCVector)cells.elementAt(AVERAGE));
    }
    public JCVector getError() {
        return new JCVector((JCVector)cells.elementAt(ERROR));
    }

    public JCVector getPressure() {
        return new JCVector((JCVector)cells.elementAt(PRESSURE));
    }

    public Vector getRow(int i) {
        Vector row = new JCVector(5);
        row.setElementAt(((JCVector)cells.elementAt(PRESSURE)).elementAt(i), i);
        row.setElementAt(((JCVector)cells.elementAt(CONFIGS)).elementAt(i), i);
        row.setElementAt(((JCVector)cells.elementAt(CURRENT)).elementAt(i), i);
        row.setElementAt(((JCVector)cells.elementAt(AVERAGE)).elementAt(i), i);
        row.setElementAt(((JCVector)cells.elementAt(ERROR)).elementAt(i), i);
        return row;
    }
    
    public void setTemp(double temp) {
        this.temp = temp;
    }
    
    public boolean addRow(int row, Object label, Vector v) {
        boolean b = super.addRow(row, label, v);
        update();
        return b;
    }

    public boolean addColumn(int col, Object label, Vector v) {
        boolean b = super.addColumn(col, label, v);
        columns = ((JCVector)cells.elementAt(PRESSURE)).size();
        update();
        return b;
    }
    
    /**
    * Sets the cell value of the JCVector "cell", abstract method from interface EditableTableData.
    */
    public boolean setTableDataItem(Object value, int row, int column) {
        //System.out.println("setTableDataItem("+row+","+column+")="+value);
        boolean result = super.setTableDataItem(getDouble(value), row, column);
        //System.out.println("NOW setTableDataItem("+row+","+column+")="+((JCVector)cells.elementAt(row)).elementAt(column)+ ((Object)((JCVector)cells.elementAt(row)).elementAt(column)).getClass());
        columns = ((JCVector)cells.elementAt(PRESSURE)).size();
        update();
        return result;
    }

    public void update() {
        this.setDataChanged(JCTblEnum.ALL, JCTblEnum.ALL, 0, 0, TableDataEvent.RESET);
    }

    /**
    * From VectorDataSource, fires a TabelDataEvent to all registered listeners by calling
    * their dataChanged method --- DualCompoundCells, DualContribVeiw and DualValueView
    * are listeners/ here.
    **/
    public void setDataChanged(int row, int column, int num_affected, int destination, int command) {
        //recalc();
        super.setDataChanged(row, column, num_affected, destination, command);
        //if (column == 1) System.out.println("setDataChanged3("+row+","+column+")="+((JCVector)cells.elementAt(row)).elementAt(column));
    }


    public Double getDouble(Object obj) {
        Double value;
        if (obj instanceof Double) {
    	value = (Double) obj;
        }
        else {
    	value = new Double(0.0);
        }
        return value;
    }

    /**
     * Gets the datasource's NumColumns value.
     * @see #setNumColumns
     */
    public int getNumColumns() {
        return columns;
    }

    /**
    * Gets the datasource's NumRows value.
     * @see #setNumRows
     */
    public int getNumRows() {
        rows = 5;
        return rows;    
    }
}
