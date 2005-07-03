package etomica.math.geometry;


/**
 * A circle in 2D.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jun 28, 2005 by kofke
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
