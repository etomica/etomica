/*
 * History
 * Created on Oct 28, 2004 by kofke
 */
package etomica.utility;

import java.lang.reflect.Array;

import etomica.Activity;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public final class Arrays {

	private Arrays() {	}

	public static Object[] resizeArray(Object[] array, int newSize) {
		Object[] newArray = (Object[])Array.newInstance(array.getClass().getComponentType(),newSize);
		int minSize = Math.min(array.length,newSize);
		System.arraycopy(array,0,newArray,0,minSize);
		return newArray;
	}
	

	public static Object[] addObject(Object[] objects, Object newObject) {
		objects = etomica.utility.Arrays.resizeArray(objects,objects.length+1);
		objects[objects.length-1] = newObject;
		return objects;
	}
	
	public static Object[] removeObject(Object[] objects, Object object) {
		int length = objects.length;
		for (int i=0; i<length; i++) {
			if (objects[i] == object) {
				Object lastObject = objects[length-1];
				objects = (Activity[])etomica.utility.Arrays.resizeArray(objects,length-1);
				if (i < length-2) {
					System.arraycopy(objects,i+1,objects,i,length-i-2);
				}
				if (i < length-1) {
					objects[length-2] = lastObject;
				}
				break;
			}
		}
		return objects;
	}

}
