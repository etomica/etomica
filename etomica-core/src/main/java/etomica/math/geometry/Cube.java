/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;
import etomica.space3d.Vector3D;

/**
 * A geometric cube.
 * 
 * @author kofke
 *  
 */
public class Cube extends Hexahedron {

    /**
     * Constructs a cube of unit size.
     */
    public Cube(Space embeddedSpace) {
        this(embeddedSpace, 1.0);
    }
    
    /**
     * Constructs a cube with edge length having the given value.
     * @param size edge length of the cube
     */
    //TODO generalize for any embedded space
    public Cube(Space embeddedSpace, double size) {
        super(embeddedSpace);
        setEdgeLength(size);
    }
    
    /**
     * Returns edgeLength^3.
     */
    public double getVolume() {
        return edgeLength*edgeLength*edgeLength;
    }
    
    /**
     * Returns 6*edgeLength^2
     */
    public double getSurfaceArea() {
        return 6*edgeLength*edgeLength;
    }
    
    /**
     * Returns 12*edgeLength
     */
    public double getPerimeter() {
        return 12*edgeLength;
    }

    /**
     * Returns <code>true</code> if the given vector lies inside (or on the surface of)
     * this cell, <code>false</code> otherwise.
     */
    public boolean contains(Vector v) {
        double x = v.getX(0)-position.getX(0);
        double y = v.getX(1)-position.getX(1);
        double z = v.getX(2)-position.getX(2);
        return (x>=n) && (x<=p) && (y>=n) && (y<=p) && (z>=n) && (z<=p);
    }
    
    /**
     * @return Returns the size, which is the length the edge of the cube.
     */
    public double getEdgeLength() {
        return edgeLength;
    }
    /**
     * @param edgeLength The size to set.
     */
    public void setEdgeLength(double edgeLength) {
        this.edgeLength = edgeLength;
        n = -0.5*edgeLength;
        p = +0.5*edgeLength;
        updateVertices();
    }
    
    public void updateVertices() {
        ((Vector3D)vertices[0]).E(n,n,n);
        ((Vector3D)vertices[1]).E(n,n,p);
        ((Vector3D)vertices[2]).E(n,p,n);
        ((Vector3D)vertices[3]).E(n,p,p);
        ((Vector3D)vertices[4]).E(p,n,n);
        ((Vector3D)vertices[5]).E(p,n,p);
        ((Vector3D)vertices[6]).E(p,p,n);
        ((Vector3D)vertices[7]).E(p,p,p);
        applyTranslationRotation();
    }
    
    private double edgeLength;
    private double n, p;//n = -size/2, p = +size/2

}
