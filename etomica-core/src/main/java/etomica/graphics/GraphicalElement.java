/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

/**
 * Interface for a simulation element that can make a graphical component
 */
public interface GraphicalElement {

    /**
        * Interface for a Simulation element that would be used in a simulation graphical user interface (GUI)
        * 
        * @param obj An object that might be used to specify the graphic that the GraphicalElement is to return.
        * In most cases the GraphicalElement ignores this parameter, and it can be set to null.
        * @return A Component that can be used in the GUI of a graphical simulation
        * @see Device
        * @see Display
        */
    public java.awt.Component graphic(Object obj);
}//end of GraphicalElement

