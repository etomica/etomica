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

/**
 * This JCChartLabel subclass takes the background color
 * of the Data item fill color that it is attached to.
 * 
 */

public class DataLabel extends JCChartLabel {

    /**
    * Set the background color of the label based on the data fill color
    * Set the data item to which the label is attached.
    * This is used when the label's AttachMethod is ATTACH_DATAINDEX.
    * @param di the data index to attach the label to
    */
    public void setDataIndex(JCDataIndex di) {
        jclass.chart.ChartDataViewSeries series = di.getSeries();
        if (series != null) {
	    setBackground(series.getStyle().getFillColor());
        }
        super.setDataIndex(di);
    }

}
