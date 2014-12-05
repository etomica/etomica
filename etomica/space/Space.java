/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IVectorMutable;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;

public abstract class Space implements java.io.Serializable, ISpace {

    protected Space() {
    }
    
    /**
     * Returns a space instance of the given dimension, which must be 1, 2, or 3.
     * 
     * @throws IllegalArgumentException if D is not 1, 2, or 3.
     */
    public static Space getInstance(int D) {
        switch(D) {
        case 1: 
            return Space1D.getInstance();
        case 2:
            return Space2D.getInstance();
        case 3:
            return Space3D.getInstance();
        default:
            throw new IllegalArgumentException("Illegal dimension for space");
        }
    }
    
    /* (non-Javadoc)
	 * @see etomica.space.ISpace#D()
	 */
    public abstract int D();

    /* (non-Javadoc)
	 * @see etomica.space.ISpace#rootD(double)
	 */
    public abstract double rootD(double a);
    
    /* (non-Javadoc)
	 * @see etomica.space.ISpace#powerD(int)
	 */
    public abstract int powerD(int a);
    
    /* (non-Javadoc)
	 * @see etomica.space.ISpace#powerD(double)
	 */
    public abstract double powerD(double a);
    
    /* (non-Javadoc)
	 * @see etomica.space.ISpace#makeVector()
	 */
    public abstract IVectorMutable makeVector();

    /* (non-Javadoc)
	 * @see etomica.space.ISpace#makeOrientation()
	 */
    public abstract IOrientation makeOrientation();
    
    /* (non-Javadoc)
	 * @see etomica.space.ISpace#makeTensor()
	 */
    public abstract Tensor makeTensor();

    /* (non-Javadoc)
	 * @see etomica.space.ISpace#makeRotationTensor()
	 */
    public abstract RotationTensor makeRotationTensor();
    
    /* (non-Javadoc)
	 * @see etomica.space.ISpace#makeArrayD(int)
	 */
    public abstract int[] makeArrayD(int i);
    /* (non-Javadoc)
	 * @see etomica.space.ISpace#makeArrayD(double)
	 */
    public abstract double[] makeArrayD(double d);

    /* (non-Javadoc)
	 * @see etomica.space.ISpace#sphereVolume(double)
	 */
    public abstract double sphereVolume(double r);

    /* (non-Javadoc)
	 * @see etomica.space.ISpace#sphereArea(double)
	 */
    public abstract double sphereArea(double r);
    

    /**
     * Returns a Vector from the space of the given dimension.
     * 
     * @throws IllegalArgumentException if D is not 1, 2, or 3.
     */
    public static IVectorMutable makeVector(int D) {
        switch(D) {
            case 1:  return new Vector1D();
            case 2:  return new Vector2D();
            case 3:  return new Vector3D();
            default: throw new IllegalArgumentException("Space.makeVector: Requested dimension not implemented");
        }
    }
    /**
     * Returns a Vector initialized to the given set of values in the array.
     * Spatial dimension of the Vector is determined by the length of a.
     * 
     * @throws IllegalArgumentException if a.length is not 1, 2, or 3.
     */
    public IVectorMutable makeVector(double[] a) {
        switch(a.length) {
            case 1:  return new etomica.space1d.Vector1D(a);
            case 2:  return new etomica.space2d.Vector2D(a);
            case 3:  return new etomica.space3d.Vector3D(a);
            default: throw new IllegalArgumentException("Space.makeVector: Requested dimension not implemented");
        }
    }
    
    /**
     * Returns a Vector initialized to the given set of values in the array (cast to double).
     * Spatial dimension of the Vector is determined by the length of a.
     */
    public IVectorMutable makeVector(int[] k) {
        double[] a = new double[k.length];
        for(int i=0; i<k.length; i++) {a[i] = k[i];}
        return makeVector(a);
    }
    
    /* (non-Javadoc)
	 * @see etomica.space.ISpace#makeVectorArray(int)
	 */
    public IVectorMutable[] makeVectorArray(int n) {
        IVectorMutable[] vectors = new IVectorMutable[n];
        for(int i=0; i<n; i++) vectors[i] = makeVector();
        return vectors;
    }
   
}//end of Space    
