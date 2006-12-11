package etomica.units;

import java.lang.reflect.*;
import java.util.LinkedList;
import java.util.Iterator;

public final class UnitFilter {

	/*
	 * This uses the filter for all dimensions from Lister.listdimensions,
	 * printing out everything to the console.
	 */

	public static void main(String[] args) {
		LinkedList dimensionlist = Lister.listdimensions();
		LinkedList unitlist = Lister.listunits();
		for (Iterator e = dimensionlist.iterator(); e.hasNext();) {
			try {
				String enext = e.next().toString();
				Class c = Class.forName(enext.toString());
				enext = enext.substring(14, enext.length());
				System.out.println("\n" + enext + ":");
				Field dfield = c.getField("DIMENSION");
				Dimension d = (Dimension) dfield.get(null);
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
			try {
				Class c = Class.forName(e.next().toString());
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

		for (Iterator e = currentunitlist.iterator(); e.hasNext();) {
			String enext = e.next().toString();
			enext = enext.substring(14, enext.length());
			System.out.println(enext);
		}
		return currentunitlist;

	}
}
