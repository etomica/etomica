/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Oct 28, 2004 by kofke, schultz
 */
package etomica.util;

import static java.util.Arrays.copyOf;

/**
 * Non-instantiable class with static utility methods for working with arrays.
 */
public final class Arrays {

	/**
	 * Private constructor to prevent instantiation.
	 */
	private Arrays() {	}

    /**
     * Returns an array formed from adding the newObject to the elements
     * in the given array.  New array is one element larger than given array,
     * and is of same type as given array.  newObject is placed at the end of 
     * the new array.
     * @param objects array with objects to be put in new array
     * @param newObject object placed at end of new array
     * @return new array
     */
	public static <T> T[] addObject(T[] objects, T newObject) {
		objects = copyOf(objects,objects.length+1);
		objects[objects.length-1] = newObject;
		return objects;
	}

	public static int[] addInt(int[] ints, int x) {
		ints = copyOf(ints, ints.length + 1);
		ints[ints.length-1] = x;
		return ints;
	}

    /**
     * Returns an array formed by removing the given object from the given array.  
     * New array is one element smaller than given array,
     * and is of same type as given array.  Order of objects in the new
     * array is unchanged from their order in the original array. 
     * If object is not in the given array,
     * method returns original array without performing any action.
     * @param array array with the objects
     * @param object object being removed
     * @return new array with the object removed
     */
	public static <T> T[] removeObject(T[] array, T object) {
		int length = array.length;
		for (int i=0; i<length; i++) {
			if (array[i] == object) {//look for object in array
				T lastObject = array[length-1];//save last object, which is about to be dropped
				array = copyOf(array,length-1);//shorten array, dropping last object
				if (i < length-2) {//overwrite target object
					System.arraycopy(array,i+1,array,i,length-i-2);
				}
				if (i < length-1) {//recover last object
					array[length-2] = lastObject;
				}
				break;
			}
		}
		return array;
	}

	public static String toString(double[] a) {
		if (a == null)
			return "null";

		int iMax = a.length - 1;
		if (iMax == -1)
			return "[]";

		StringBuilder b = new StringBuilder();
		b.append('[');
		for (int i = 0; ; i++) {
			b.append(String.valueOf(a[i]));
			if (i == iMax)
				return b.append(']').toString();
			b.append(", ");
		}
	}

}
