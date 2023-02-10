/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

public class ParticleRecord {
    /** The index of the block that the particle is within. */
    int ijk;
    /** The number of particle within its block. */
    int l;
    /** The x-index of the block. */
    int di;
    /** The y-index of the block. */
    int dj;
    /** The z-index of the block. */
    int dk;
}
