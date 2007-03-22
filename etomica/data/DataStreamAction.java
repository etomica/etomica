package etomica.data;

import etomica.action.Action;

public abstract class DataStreamAction implements Action, java.io.Serializable {

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
        else if (obj instanceof DataProcessor) {
            actionPerformed(((DataProcessor)obj).getDataSink());
        }
    }
    
    public abstract void dataWalkerAction(Object obj);
    
    private Object start;
}
