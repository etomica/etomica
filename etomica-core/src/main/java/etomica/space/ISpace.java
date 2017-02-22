/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IVectorMutable;

public interface ISpace {

	/**
	 * The dimension of this space.
	 */
	public abstract int D();

	/**
	 * Returns the given value raised to the power 1/D, where D is the dimension of the space.
	 */
	public abstract double rootD(double a);

	/**
	 * Returns the given value raised to the Dth power, where D is the dimension of the space.
	 */
	public abstract int powerD(int a);

	/**
	 * Returns the given value raised to the Dth power, where D is the dimension of the space.
	 */
	public abstract double powerD(double a);

	/**
	 * Constructs and returns a new Vector appropriate to the space.
	 */
	public abstract IVectorMutable makeVector();

	public abstract IVectorMutable makeVector(double[] a);

	public abstract IVectorMutable makeVector(int[] k);

	/**
	 * Constructs and returns a new Orientation appropriate to the space.
	 */
	public abstract IOrientation makeOrientation();

	/**
	 * Constructs and returns a new Tensor appropriate to the space.
	 */
	public abstract Tensor makeTensor();

	/**
	 * Constructs and returns a new RotationTensor appropriate to the space.
	 */
	public abstract RotationTensor makeRotationTensor();

	/**
	 * Returns an array of dimension D, with each element equal to the given value.
	 */
	public abstract int[] makeArrayD(int i);

	/**
	 * Returns an array of dimension D, with each element equal to the given value.
	 */
	public abstract double[] makeArrayD(double d);

	/**
	 * Returns the "volume" of the "sphere" defined in the D-dimensional space.
	 * In 1-D, this is twice the radius; in 2-D the area of the circle; 
	 * in 3-D the volume of the sphere.
	 *
	 * @param r the radius
	 * @return the volume
	 */
	public abstract double sphereVolume(double r);

	/**
	 * Returns the surface "area" of the "sphere" defined in the D-dimensional space.
	 * In 1-D this is zero; in 2-D the circumference of the circle; in 3-D 
	 * the surface area of the sphere.
	 *
	 * @param r the radius
	 * @return the area
	 */
	public abstract double sphereArea(double r);

	/**
	 * Instance methods that makes and returns an array of vectors having the
	 * given number of elements.
	 * @param n number of vectors in the returned array
	 * @return an array of n new vectors made by the space instance
	 */
	public abstract IVectorMutable[] makeVectorArray(int n);
	
}