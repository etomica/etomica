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

import java.awt.Event;
import java.util.Vector;
import jclass.base.Border;
import jclass.chart.EventTrigger;
import jclass.chart.JCChart;
import jclass.chart.JCChartLabel;
import jclass.chart.JCDataIndex;
import jclass.table3.TableDataEvent;
import jclass.table3.TableDataListener;
import jclass.util.JCVector;
import jclass.chart.JCLabelGenerator;

/**
 * This class is notified when a CompoundModel instance
 * changes and updates the Chart Label text to reflect
 * a pretty version of the last value in the CompoundModel.
 * The value array represents the accumulated value of
 * each client's contributions, so this label represents
 * the final result of allowing the account to accumulate.
 * The class does double duty by formating chart JCAxis
 * annotations by implementing the JCLabelGenerator interface
 * makeLabe().
 */

public class LabelLast extends DataLabel 
                   implements TableDataListener, JCLabelGenerator {

    double y[] = {500};


    public LabelLast() {
        JCVector pressure = new JCVector();
        JCVector configs = new JCVector();
        JCVector current = new JCVector();
        JCVector average = new JCVector();
        JCVector error = new JCVector();
        for (int i=0; i<y.length; i++) {
	        pressure.addElement(new Double(y[i]));
	        configs.addElement(new Double(i));
	        current.addElement(new Double(i));
	        average.addElement(new Double(i));
	        error.addElement(new Double(i));
        }
        EnergyDataModel model = new EnergyDataModel(pressure, configs, current, average, error);
        addListener(model);
    }

    public LabelLast(EnergyDataModel one) {
        addListener(one);
    }

    public void addListener(EnergyDataModel model) {
        model.addTableDataListener(this);
    }

    /**
    * When the model data changes, dataChanged() will
    * be called. Just let our observers(the chart)
    * know so that they can re-extract the new values.
    */
    public void dataChanged(TableDataEvent ev) {
        //System.out.println("LabelLast.dataChanged "+ev.getSource());
        Object source = ev.getSource();
        if (source instanceof EnergyDataModel) {
    		JCVector average = ((EnergyDataModel)source).getAverage();
    		Double last = (Double) average.getLast();
	    	//setText("$"+formatLabel(last.doubleValue(), 0));
		    setText(formatLabel(last.doubleValue(), 0));
        }
    }

    /**
    * Format the value as an easy-to-scan way
    * string by grouping 000s with a ','
    * End result: turn "921217.0371352091" to "921,217.04"
    */

    public static String formatLabel(double value, int decimals) {
        String nice = new String();
        java.text.NumberFormat formatter = 
        java.text.NumberFormat.getCurrencyInstance();
    
        formatter.setMaximumFractionDigits( decimals );
        nice = formatter.format( value );
        /*
        String tri;
        int simple = (int) value;
        if (simple == 0) return new String("0");
        int rem, maj;
        boolean first = true;
        for (int i = 1000000000; i>999; i*=0.001) {
    		if (first && i > simple) continue;
    		rem = simple % i;
    		maj = simple - rem;
    		simple -= maj;
    		maj = (int) (maj / i);
    		tri = "" + maj;
    		if (!first) {
    			for (int j=0; j<3-tri.length(); j++) nice = nice + "0";
    		}
    		nice = nice + tri + ",";
    		//nice = nice + maj + ",";
    		first = false;
        }
        tri = "" + (simple % 1000);
        if (!first) for (int j=0; j<3-tri.length(); j++) nice = nice + "0";
        nice = nice + tri;
        if (decimals > 0) {
    		// Untested...
    		double fract = value - (double) ((int)value);
    		String small = jclass.chart.JCChartUtil.format(fract, 0);
    		// Not i18n! Would need to be much more careful if this wasn't
    		// just a demo...
    		nice = nice + small.substring(small.indexOf("."));
        }
        */
        return nice;
    }
    
    /**
    * Callback method called by JCChart when a label is needed. 
    * Create a ChartText instance representing the supplied
    * value multiplied by the scale factor.
    * @param value axis value that is to be labelled.
    * @param precision numeric precision to be used.
    */
    public Object makeLabel(double value, int precision) {
        jclass.chart.ChartText c = new jclass.chart.ChartText();
        String label = formatLabel(value, 0);
        if (label.equals("0")) { 
    		label = "<" + label + ">"; 
    	}
        c.setText(label);   
    	//System.out.println("formatted:("+value+","+precision+") to "+c.getText());
        return label;
    }

}
