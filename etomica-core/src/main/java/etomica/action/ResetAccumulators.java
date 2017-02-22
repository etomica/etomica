/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
