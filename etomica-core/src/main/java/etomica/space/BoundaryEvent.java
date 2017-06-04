/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IBoundaryEvent;

public class BoundaryEvent implements IBoundaryEvent {

    protected Boundary boundary = null;
    
    public BoundaryEvent(Boundary _boundary) {
        boundary = _boundary;
    }
    
    public Boundary getBoundary() {
        return boundary;
    }
}
