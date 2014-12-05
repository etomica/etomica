/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.integrator.IntegratorHard;

/**
 * Interface for a meter that can return a value based on a hard-potential collision.
 */
public interface DataSourceCollisional extends IntegratorHard.CollisionListener {
    public double collisionValue(IntegratorHard.Agent agent);
}