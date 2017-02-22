/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Color;
import java.awt.TextArea;

import javax.swing.JFrame;

import etomica.action.IAction;
import etomica.data.AccumulatorHistory;
import etomica.data.IData;
import etomica.data.types.DataFunction.DataInfoFunction;

/**
 * Action that opens a new window and dumps the coordinates into the window.
 * @author Andrew Schultz
 */
public class ActionHistoryWindow implements IAction {
    protected final AccumulatorHistory accumulatorHistory;
    
    public ActionHistoryWindow(AccumulatorHistory accumulatorHistory) {
        this.accumulatorHistory = accumulatorHistory;
    }
    
    public void actionPerformed() {
        JFrame f = new JFrame();
        TextArea textArea = new TextArea();
        textArea.setEditable(false);
        textArea.setBackground(Color.white);
        textArea.setForeground(Color.black);
        IData data = accumulatorHistory.getData();
        IData xData = ((DataInfoFunction)accumulatorHistory.getDataInfo()).getXDataSource().getIndependentData(0);
        for (int i=0; i<data.getLength(); i++) {
            double x = xData.getValue(i);
            double y = data.getValue(i);
            if (Double.isNaN(x) || Double.isNaN(y)) continue;
            textArea.append(x+"\t"+y+"\n");
        }
        f.add(textArea);
        f.pack();
        f.setSize(400,600);
        f.setVisible(true);
    }
}