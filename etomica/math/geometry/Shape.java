package etomica.math.geometry;

import etomica.space.IVector;

/**
 * Interface for an object that defines a closed region of space.
 * 
 * @author David Kofke
 *  
 */
public interface Shape {

    /**
     * Returns the (hyper-)volume of the enclosed region. In 1D this is a
     * length, in 2D it is an area, in 3D it is a volume.
     */
    public double getVolume();

    /**
     * Returns true if the point described by the vector lies in or on the
     * enclosed region.
     */
    public boolean contains(IVector r);
    
    /**
     * Sets the position of the shape, which typically represents some central
     * point within it.  Specific definition depends on the definition of the shape.
     */
    public void setPosition(IVector r);

    /**
     * Returns the spatial dimension of the shape. For example, D = 2 for a
     * square, D = 3 for a cube.
     */
    public int D();
}
