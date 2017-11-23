/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.action.IAction;

public abstract class DataStreamAction implements IAction, java.io.Serializable {

    public DataStreamAction() {
        super();
    }

    public void setStart(Object start) {
        this.start = start;
    }
    
    public void actionPerformed() {
        actionPerformed(start);
    }
    
    public void actionPerformed(Object obj) {
        dataWalkerAction(obj);
        if (obj instanceof DataPipeForked) {
            IDataSink[] dataSinks = ((DataPipeForked)obj).getDataSinks();
            for (int i=0; i<dataSinks.length; i++) {
                actionPerformed(dataSinks[i]);
            }
            if (obj instanceof AccumulatorAverageFixed) {
                actionPerformed(((AccumulatorAverageFixed) obj).getBlockDataSink());
            }
        }
        else if (obj instanceof DataSplitter) {
            int n = ((DataSplitter)obj).getNumDataSinks();
            for (int i=0; i<n; i++) {
                actionPerformed(((DataSplitter)obj).getDataSink(i));
            }
        }
        else if (obj instanceof DataGroupSplitter) {
            int n = ((DataGroupSplitter)obj).getNumDataSinks();
            for (int i=0; i<n; i++) {
                actionPerformed(((DataGroupSplitter)obj).getDataSink(i));
            }
        }
        else if (obj instanceof DataProcessor) {
            actionPerformed(((DataProcessor)obj).getDataSink());
        }
    }
    
    public abstract void dataWalkerAction(Object obj);
    
    private Object start;
}
