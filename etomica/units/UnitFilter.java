package etomica.units;

import java.lang.reflect.*;
import java.util.LinkedList;
import java.util.Iterator;

public final class UnitFilter {

	/*
	 * This uses the filter for all dimensions from Lister.listdimensions,
	 * printing out everything to the console.
	 */

	public static final void main(String[] args) {
		LinkedList<String> dimensionlist = Lister.listdimensions();
		for (Iterator e = dimensionlist.iterator(); e.hasNext();) {
			try {
				String enext = e.next().toString();
				System.out.println("\n" + enext + ":");
				Class c = Class.forName("etomica.units." + enext);
				Field dfield = c.getField("DIMENSION");
				Dimension d = (Dimension) dfield.get(null);
				filter(d);
			} catch (Throwable er) {
				System.err.println(er);
			}
		}
	}

	/*
	 * This filter accepts a Dimension, goes through all units obtained from
	 * Lister.listunits, then removes any that are not of the given dimension.
	 */

	public static LinkedList<String> filter(Dimension dim) {

		LinkedList<String> unitlist = Lister.listunits();

		for (Iterator e = unitlist.iterator(); e.hasNext();) {
			try {
				Class c = Class.forName("etomica.units." + e.next());
				Field ufield = c.getField("UNIT");
				Unit currentunit = (Unit) ufield.get(null);
				if (currentunit.dimension() != dim) {
					e.remove();
				}
			} catch (Throwable er) {
				System.err.println(er);
				System.out.println(e.next());
			}
		}

		for (Iterator e = unitlist.iterator(); e.hasNext();) {
			System.out.println(e.next());
		}
		return unitlist;

	}
}
