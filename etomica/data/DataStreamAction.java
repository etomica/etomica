package etomica.data;

import etomica.api.IAction;

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
            DataSink[] dataSinks = ((DataPipeForked)obj).getDataSinks();
            for (int i=0; i<dataSinks.length; i++) {
                actionPerformed(dataSinks[i]);
            }
        }
        else if (obj instanceof DataSplitter) {
            int n = ((DataSplitter)obj).getNumDataSinks();
            for (int i=0; i<n; i++) {
                actionPerformed(((DataSplitter)obj).getDataSink(i));
            }
        }
        else if (obj instanceof DataProcessor) {
            actionPerformed(((DataProcessor)obj).getDataSink());
        }
    }
    
    public abstract void dataWalkerAction(Object obj);
    
    private Object start;
}
