package etomica.action;

import etomica.data.AccumulatorHistory;
import etomica.util.HistoryScrolling;

public class ResetAccumulatorsAveraged extends ResetAccumulators {

    public void dataWalkerAction(Object obj) {
        if (obj instanceof AccumulatorHistory &&
            ((AccumulatorHistory)obj).getHistory() instanceof HistoryScrolling) {
            return;
        }
        super.dataWalkerAction(obj);
    }
    private static final long serialVersionUID = 1L;

}
