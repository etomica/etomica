package etomica.action;

import etomica.data.DataAccumulator;
import etomica.data.DataStreamAction;

public class ResetAccumulators extends DataStreamAction {

    public ResetAccumulators() {
        super();
    }

    public void dataWalkerAction(Object obj) {
        if (obj instanceof DataAccumulator) {
            ((DataAccumulator)obj).reset();
        }
    }

}
