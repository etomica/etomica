/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IBoundary;
import etomica.api.IBoundaryEvent;

public class BoundaryEvent implements IBoundaryEvent {

    protected IBoundary boundary = null;
    
    public BoundaryEvent(IBoundary _boundary) {
        boundary = _boundary;
    }
    
    public IBoundary getBoundary() {
        return boundary;
    }
}
