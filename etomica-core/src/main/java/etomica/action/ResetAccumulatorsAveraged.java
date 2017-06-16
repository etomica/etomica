/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.data.AccumulatorHistory;
import etomica.data.history.HistoryScrolling;

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
