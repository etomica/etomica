/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import java.io.File;
import java.io.FileInputStream;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.StringTokenizer;
import java.util.jar.JarEntry;
import java.util.jar.JarInputStream;

import etomica.units.dimensions.Dimension;
import etomica.units.systems.UnitSystem;

public final class Lister {

	/*
	 * This method takes the classes form the unitclasses function and removes
	 * any that don't have the Dimension superclass, Null, Undefined, and any
	 * dimensions with "dimension" in the title.
	 */

	public static LinkedList<String> listdimensions() {
		LinkedList<String> dimensionlist = classlist();
		for (Iterator<String> e = dimensionlist.iterator(); e.hasNext();) {
			try {
				Class c = Class.forName(e.next().toString());
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

	public static LinkedList<String> listUnitSystems() {
		LinkedList<String> unitsystemlist = classlist();
		for (Iterator<String> e = unitsystemlist.iterator(); e.hasNext();) {
			try {
				Class c = Class.forName(e.next().toString());
				if (c.getSuperclass() != UnitSystem.class) {
					e.remove();
				}
			} catch (ClassNotFoundException ee) {
				System.out.println(ee);
			} catch (Exception eee) {
				System.out.println(eee);
			}
		}

		// for (int i = 0; i < unitsystemlist.size(); i++) {
		// System.out.println(unitsystemlist.get(i).toString());
		// }
		// System.out.println();
		return unitsystemlist;
	}
	
	/*
	 * This method takes the classes form the unitclasses function and removes
	 * any that don't have a SimpleUnit or CompoundUnit superclass, and removes
	 * the pixel class and any containing the word unit.
	 */

	public static LinkedList<String> listUnits() {
		LinkedList<String> unitslist = classlist();
		for (Iterator<String> e = unitslist.iterator(); e.hasNext();) {
			try {
				Class c = Class.forName(e.next().toString());
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

	public static LinkedList<String> listUnitsInSystem(String unitSystem){
		LinkedList<String> unitList = listUnits();
		LinkedList<String> unitsInSystem = new LinkedList<String>();
		
		try {
			UnitSystem us = UnitFilter.stringToUnitSystem(unitSystem);
			Class c = us.getClass();
			Method[] methodSet = c.getDeclaredMethods();
			for (int i=0; i<methodSet.length;i++){
				String s = methodSet[i].invoke(us).getClass().getName();
				//System.out.println(s);
				if (unitList.contains(s)){
					unitsInSystem.add(s);
				}
			}
		} catch (Throwable er) {
			//System.err.println(er);
		}
		/*for(Iterator e = unitsInSystem.iterator(); e.hasNext();){
			System.out.println(e.next().toString());
		}*/
		return unitsInSystem;
	}
	public static String[] uisString(String us) {
		LinkedList<String> dl = listUnitsInSystem(us);
		String[] ds = new String[dl.size()];
		int x = 0;
		for (Iterator<String> e = dl.iterator(); e.hasNext();) {
			try {
				String enext = e.next().toString();
				enext = enext.substring(14, enext.length());
				//System.out.println(enext);
				ds[x] = enext;
				x += 1;
			} catch (Throwable er) {
				System.err.println(er);
			}

		}
		return ds;
	}
	
	public static String[] dimensionString() {
		LinkedList<String> dl = listdimensions();
		String[] ds = new String[dl.size()];
		int x = 0;
		for (Iterator<String> e = dl.iterator(); e.hasNext();) {
			try {
				String enext = e.next().toString();
				enext = enext.substring(14, enext.length());
				ds[x] = enext;
				x += 1;
			} catch (Throwable er) {
				System.err.println(er);
			}

		}
		return ds;
	}

	public static String[] unitSystemString() {
		LinkedList<String> ul = listUnitSystems();
		String[] us = new String[ul.size()];
		int x = 0;
		for (Iterator<String> e = ul.iterator(); e.hasNext();) {
			try {
				String enext = e.next().toString();
				enext = enext.substring(22, enext.length());
				us[x] = enext;
				x += 1;
			} catch (Throwable er) {
				System.err.println(er);
			}

		}
		return us;
	}

	public static String[] unitString() {
		LinkedList<String> ul = listUnits();
		String[] us = new String[ul.size()];
		int x = 0;
		for (Iterator<String> e = ul.iterator(); e.hasNext();) {
			try {
				String enext = e.next().toString();
				enext = enext.substring(14, enext.length());
				us[x] = enext;
				x++;
			} catch (Throwable er) {
				System.err.println(er);
			}

		}
		return us;
	}

	
	/*
	 * This method finds all the files in the src/etomica/units directory, takes
	 * any with a .java extension, removes the .java, and adds them all to a
	 * linked list of strings.
	 * 
	 * Replaced with classlist();
	 */

/*	public static LinkedList unitsclasses() {
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
				clist.add(currentitem); // System.out.println(currentitem);
			}
		}
		return clist;
	}*/

	/*
	 * This replaces unitclasses.
	 */
	public static LinkedList<String> classlist() {

		LinkedList<String> clist = new LinkedList<String>();
		ArrayList<String> paths = new ArrayList<String>();

		String cpath = System.getProperty("java.class.path");
		// System.out.println( "From Etomica plugin: java.class.path = " + cpath
		// );
		StringTokenizer tz = new StringTokenizer(cpath, PATHSEP);
		while (tz.hasMoreTokens()) {
			paths.add(tz.nextToken());
		}
		for (int i = 0; i < paths.size(); i++) {
			File file = new File(paths.get(i));
			if (file.isDirectory()) {
				if (!processDirectory(file, file.getAbsolutePath()).isEmpty())
					clist.addAll(processDirectory(file, file
									.getAbsolutePath()));
			} else {
				if (!processJar(file, file.getAbsolutePath()).isEmpty())
					clist.addAll(processJar(file, file.getAbsolutePath()));
			}

		}

		return clist;
	}

	/**
	 * Recursively store all class names found in this directory and its
	 * subdirectories
	 */
	private static LinkedList<String> processDirectory(File directory, String basepath) {
		// System.out.println("Processing " + directory.getPath());
		LinkedList<String> addlist = new LinkedList<String>();
		File[] files = directory.listFiles();
		for (int i = 0; i < files.length; i++) {
			File currentFile = files[i];
			String name = currentFile.getName();

			if (currentFile.isDirectory()) {
				if (!processDirectory(currentFile, basepath).isEmpty())
					addlist.addAll(processDirectory(currentFile, basepath));

			} else if (name.endsWith(".class")) {
				if (!add(name, currentFile.getPath(), basepath).isEmpty())
					addlist.addAll(add(name, currentFile.getPath(), basepath));

			} else if (name.endsWith(".jar")) {
				if (!processJar(currentFile, basepath).isEmpty())
					addlist.addAll(processJar(currentFile, basepath));
			}
		}
		return addlist;
	}

	/**
	 * Store all class names found in this jar file
	 */
	private static LinkedList<String> processJar(File jarfile, String basepath) {
		// System.out.println("Processing JAR " + jarfile.getPath());
		LinkedList<String> addlist = new LinkedList<String>();
		try {
			JarInputStream jarIS = new JarInputStream(new FileInputStream(
					jarfile));

			JarEntry entry = null;
			while ((entry = jarIS.getNextJarEntry()) != null) {
				String name = entry.getName();
				if (name.endsWith(".class")) {
					// System.out.println( entry.getAttributes().toString() );
					if (!add(name, jarfile.getPath(), basepath).isEmpty())
						addlist.addAll(add(name, jarfile.getPath(), basepath));
				}
			}
		} catch (Exception e) {
		}
		return addlist;
	}

	private static LinkedList<String> add(String origname, String path, String basepath) {

		LinkedList<String> addlist = new LinkedList<String>();
		String name = origname;
		if (path.startsWith(basepath)) {
			try {
				name = path.substring(basepath.length() + 1);
			} catch (Exception e) {
			}
		}

		if (!name.startsWith("etomica") || name.indexOf('$') >= 0)
			return addlist;

		// Get class object
		name = name.replace('\\', '.');
		name = name.replace('/', '.');

		if (name.startsWith("etomica.units")) {

			if (name.endsWith(".class")) {
				name = name.substring(0, name.length() - 6);
			}
			Class thisclass = null;
			// Pressure2D and Area don't have units yet.
			//if (name.indexOf("Pressure2D") > -1
			//		|| name.toString().indexOf("Area") > -1) {
			//} else {
				try {
					// System.err.println( "Trying class " + name );
					thisclass = Class.forName(name);
					addlist.add(name);
				} catch (java.lang.ClassFormatError e) {
					System.out.println("Could not access class " + name + ": "
							+ e.getLocalizedMessage());
					return addlist;
				} catch (java.lang.VerifyError e) {
					System.out.println("Could not access class " + name + ": "
							+ e.getLocalizedMessage());
					return addlist;
				} catch (java.lang.NoClassDefFoundError e) {
					System.out.println("Could not access class " + name + ": "
							+ e.getLocalizedMessage());
					return addlist;
				} catch (Exception e) {
					System.out.println("Could not access class " + name + ": "
							+ e.getLocalizedMessage());
					return addlist;
				} catch (ExceptionInInitializerError e) {
					System.out.println("Could not access class " + name + ": "
							+ e.getLocalizedMessage());
					return addlist;
				}
				if ((thisclass.getModifiers() & java.lang.reflect.Modifier.ABSTRACT) != 0) {
					return addlist;
				}
			}
		//}
		return addlist;
	}

	private static String PATHSEP = System.getProperty("path.separator");
}
