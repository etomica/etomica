/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

public class MCMoveTrialCompletedEvent extends MCMoveEvent implements java.io.Serializable {

    public MCMoveTrialCompletedEvent(MCMoveManager moveManager, boolean accepted) {
        super();
        isAccepted = accepted;
        this.moveManager = moveManager;
    }
    
    public MCMove getMCMove() {
        return moveManager.getSelectedMove();
    }
    
    public boolean isAccepted() {
        return isAccepted;
    }
    
    private final boolean isAccepted;
    private final MCMoveManager moveManager;
    public double chi;
}
