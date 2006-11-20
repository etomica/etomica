package etomica.units;

import java.io.*;
import java.util.LinkedList;
import java.util.Iterator;

public final class Lister {

	/*
	 * This method takes the classes form the unitclasses function and removes
	 * any that don't have the Dimension superclass, Null, Undefined, and any
	 * dimensions with "dimension" in the title.
	 */

	public static LinkedList listdimensions() {
		LinkedList dimensionlist = unitsclasses();
		for (Iterator e = dimensionlist.iterator(); e.hasNext();) {
			try {
				Class c = Class.forName("etomica.units." + e.next());
				if (c.getSuperclass() != Dimension.class
						|| c.toString().indexOf("Dimension") > -1
						|| c.toString().indexOf("Null") > -1
						|| c.toString().indexOf("Undefined") > -1) {
					e.remove();
				}
			} catch (ClassNotFoundException ee) {
				System.out.println(ee);
			} catch (Exception eee) {
				System.out.println(eee);
			}
		}
		// for (int i = 0; i < dimensionlist.size(); i++) {
		// System.out.println(dimensionlist.get(i).toString());
		// }
		// System.out.println();
		return dimensionlist;
	}

	/*
	 * This method takes the classes form the unitclasses function and removes
	 * any that don't have a SimpleUnit or CompoundUnit superclass, and removes
	 * the pixel class and any containing the word unit.
	 */

	public static LinkedList listunits() {
		LinkedList unitslist = unitsclasses();
		for (Iterator e = unitslist.iterator(); e.hasNext();) {
			try {
				Class c = Class.forName("etomica.units." + e.next());
				if ((c.getSuperclass() != SimpleUnit.class && c.getSuperclass() != CompoundUnit.class)
						|| c.toString().indexOf("Unit") > -1
						|| c.toString().indexOf("Pixel") > -1) {
					e.remove();
				}
			} catch (ClassNotFoundException ee) {
				System.out.println(ee);
			} catch (Exception eee) {
				System.out.println(eee);
			}
		}
		/*
		 * for (int i = 0; i < unitslist.size(); i++) {
		 * System.out.print(unitslist.get(i).toString() + "\n"); }
		 * System.out.println();
		 */
		return unitslist;
	}

	/*
	 * This method finds all the files in the src/etomica/units directory, takes
	 * any with a .java extension, removes the .java, and adds them all to a
	 * linked list of strings.
	 */
	public static LinkedList unitsclasses() {
		File dir = new File("src/etomica/units");
		String[] listy = dir.list();
		LinkedList clist = new LinkedList();
		for (int i = 0; i < listy.length; i++) {
			String currentitem = listy[i];
			int dotspot = currentitem.indexOf(".");
			String filetype = currentitem.subSequence(dotspot + 1,
					currentitem.length()).toString();
			if (filetype.equals("java")) {
				currentitem = currentitem.substring(0, dotspot);
				clist.add(currentitem);
				// System.out.println(currentitem);
			}
		}
		return clist;
	}

}
