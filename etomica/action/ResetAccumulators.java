package etomica.action;

import etomica.data.DataAccumulator;
import etomica.data.DataStreamAction;

public class ResetAccumulators extends DataStreamAction {

    public void dataWalkerAction(Object obj) {
        if (obj instanceof DataAccumulator) {
            ((DataAccumulator)obj).reset();
        }
    }

    private static final long serialVersionUID = 1L;
}
