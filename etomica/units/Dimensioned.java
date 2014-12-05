/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

/**
 * Interface for an object with an associated physical dimension.
 * Indicates that the object has methods for setting and getting its units.
 */

public interface Dimensioned {
    
    public void setUnit(Unit u);
    public Unit getUnit();
    public Dimension getDimension();
    
}