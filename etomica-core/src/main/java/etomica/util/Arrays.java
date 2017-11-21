/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Oct 28, 2004 by kofke, schultz
 */
package etomica.util;

import java.util.LinkedList;

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
	public static Object[] addObject(Object[] objects, Object newObject) {
		objects = copyOf(objects,objects.length+1);
		objects[objects.length-1] = newObject;
		return objects;
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
	public static Object[] removeObject(Object[] array, Object object) {
		int length = array.length;
		for (int i=0; i<length; i++) {
			if (array[i] == object) {//look for object in array
				Object lastObject = array[length-1];//save last object, which is about to be dropped
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

    /**
     * Converts the given array of integers into a string with a comma-separated
     * list of the integers enclosed by braces.  Null argument or empty array
     * are converted to the string "{}".
     */
    public static String toString(int[] array) {
        if(array == null) return "{}";
        StringBuffer string = new StringBuffer("{");
        for(int i=0; i<array.length-1; i++) {
            string.append(Integer.toString(array[i])).append(",");
        }
        if(array.length > 0) string.append(Integer.toString(array[array.length-1]));
        string.append("}");
        return string.toString();
    }

    /**
     * Converts the given array of doubles into a string with a comma-separated
     * list of the doubles enclosed by braces.  Null argument or empty array
     * are converted to the string "{}".
     */
    public static String toString(double[] array) {
        if(array == null) return "{}";
        StringBuffer string = new StringBuffer("{");
        for(int i=0; i<array.length-1; i++) {
            string.append(Double.toString(array[i])).append(",");
        }
        if(array.length > 0) string.append(Double.toString(array[array.length-1]));
        string.append("}");
        return string.toString();
    }
    /**
     * Converts the given array of floats into a string with a comma-separated
     * list of the floats enclosed by braces.  Null argument or empty array
     * are converted to the string "{}".
     */
    public static String toString(float[] array) {
        if(array == null) return "{}";
        StringBuffer string = new StringBuffer("{");
        for(int i=0; i<array.length-1; i++) {
            string.append(Float.toString(array[i])).append(",");
        }
        if(array.length > 0) string.append(Float.toString(array[array.length-1]));
        string.append("}");
        return string.toString();
    }
    /**
     * Converts the given array of booleans into a string with a comma-separated
     * list of the booleans enclosed by braces.  Null argument or empty array
     * are converted to the string "{}".
     */
    public static String toString(boolean[] array) {
        if(array == null) return "{}";
        StringBuffer string = new StringBuffer("{");
        for(int i=0; i<array.length-1; i++) {
            string.append(Boolean.toString(array[i])).append(",");
        }
        if(array.length > 0) string.append(Boolean.toString(array[array.length-1]));
        string.append("}");
        return string.toString();
    }

    /**
     * Creates a new LinkedList filled with the elements of the given array, in the same order.
     */
    public static LinkedList toList(Object[] obj) {
        LinkedList list = new LinkedList();
        for(int i=0; i<obj.length; i++) {
            list.add(obj[i]);
        }
        return list;
    }
}
