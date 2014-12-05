/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

public class MCMoveTrialInitiatedEvent extends MCMoveEvent implements java.io.Serializable {

    public MCMoveTrialInitiatedEvent(MCMoveManager moveManager) {
        super();
        this.moveManager = moveManager;
    }
    
    public MCMove getMCMove() {
        return moveManager.getSelectedMove();
    }
    
    private final MCMoveManager moveManager;
}
