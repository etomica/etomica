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
            int n = ((DataPipeForked)obj).getDataSinkCount();
            for (int i=0; i<n; i++) {
                actionPerformed(((DataPipeForked)obj).getDataSink(i));
            }
        }
        else if (obj instanceof DataProcessor) {
            actionPerformed(((DataProcessor)obj).getDataSink());
        }
    }
    
    public abstract void dataWalkerAction(Object obj);
    
    public String getLabel() {
        return label;
    }
    
    public void setLabel(String label) {
        this.label = label;
    }
    
    private Object start;
    private String label;
}
