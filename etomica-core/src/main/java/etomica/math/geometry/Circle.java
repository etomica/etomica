/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;


/*
 * History Created on Jun 28, 2005 by kofke
 */

/**
 * A circle in 2D.
 * 
 * @author David Kofke
 *
 */
public class Circle extends Hypersphere {

    /**
     * Creates circle with unit radius
     */
    public Circle() {
        this(1.0);
    }

    /**
     * Creates circle of the given radius.
     *
     * @param radius the radius of the circle
     */
    public Circle(double radius) {
        super(2, radius);
    }
    
    /**
     * Returns the volume for the present sphere radius.
     */
    public double getVolume() {
        return Math.PI * radius * radius;
    }

}
