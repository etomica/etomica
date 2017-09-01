/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.gui;

import java.lang.reflect.Field;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.units.Unit;
import etomica.units.dimensions.Dimension;
import etomica.units.systems.UnitSystem;

public final class UnitFilter {

	/*
	 * This uses the filter for all dimensions from Lister.listdimensions,
	 * printing out everything to the console.
	 */

	public static void main(String[] args) {
		LinkedList dimensionlist = Lister.listdimensions();
		LinkedList unitlist = Lister.listUnits();
		for (Iterator e = dimensionlist.iterator(); e.hasNext();) {
			try {
				String enext = e.next().toString();
				Dimension d = stringToDim(enext);
				enext = enext.substring(14, enext.length());
				System.out.println("\n" + enext + ":");
				filter(d, unitlist);
			} catch (Throwable er) {
				System.err.println(er);
			}
		}

	}

	/*
	 * This filter accepts a Dimension and a list to filter, goes through all
	 * units obtained from Lister.listunits, then removes any that are not of
	 * the given dimension.
	 */

	public static LinkedList filter(Dimension dim, LinkedList unitlist) {
		LinkedList currentunitlist = (LinkedList) unitlist.clone();
		for (Iterator e = currentunitlist.iterator(); e.hasNext();) {
			try {String enext = e.next().toString();
				enext = checkString(enext);
				Class c = Class.forName(enext);
				Field ufield = c.getField("UNIT");
				Unit currentunit = (Unit) ufield.get(null);
				if (currentunit.dimension() != dim) {
					e.remove();
				}
			} catch (Throwable er) {
//				System.err.println(er);
//				System.out.println(e.next());
			}
		}

		for (Iterator e = currentunitlist.iterator(); e.hasNext();) {
			String enext = e.next().toString();
			enext = checkString(enext);
			enext = enext.substring(14, enext.length());
//			System.out.println(enext);
		}
		return currentunitlist;

	}
	public static String[] filter(Dimension dim, String[] unitString){
		LinkedList ll = new LinkedList();
		String[] newUnitString;
		for (int i = 0; i<unitString.length;i++){
		ll.add(unitString[i]);
		}
		ll = filter(dim,ll);
		newUnitString = new String[ll.size()];
		int count = 0;
		for (Iterator e = ll.iterator(); e.hasNext();){
			String enext = e.next().toString();
			newUnitString[count]=enext;
			count++;
			}
		return newUnitString;
	}
	public static String checkString(String string){
		if (string.length()>13){
			if (!string.substring(0,13).equals("etomica.units.")){
				string = "etomica.units." + string;}}
			else {
				string = "etomica.units." + string;
			}
		return string;
	}
	public static String checkUnitSystemString(String string){
		if (string.length()>20){
			if (!string.substring(0,20).equals("etomica.units.systems.")){
				string = "etomica.units.systems." + string;}}
			else {
				string = "etomica.units.systems." + string;
			}
		return string;
	}
	public static Dimension stringToDim(String string){
		String dimString = checkString(string);
		
		try {
			Class c = Class.forName(dimString);
			Field dfield = c.getField("DIMENSION");
			Dimension d = (Dimension) dfield.get(null);
			return d;
		} catch (Throwable er) {
			//System.err.println(er);
		}
		return null;
	}
	public static UnitSystem stringToUnitSystem(String string){
		String unitSystemString = checkUnitSystemString(string);
		
		try {
			Class c = Class.forName(unitSystemString);
			
			UnitSystem u = (UnitSystem)c.newInstance();
			return u;
		} catch (Throwable er) {
			//System.err.println(er);
		}
		return null;
	}	public static Unit stringToUnit(String string){
		String unitString = checkString(string);
		
		try {
			Class c = Class.forName(unitString);
			Field ufield = c.getField("UNIT");
			Unit u = (Unit) ufield.get(null);
			return u;
		} catch (Throwable er) {
			//System.err.println(er);
		}
		return null;
	}
}
