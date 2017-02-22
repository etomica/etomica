The sources here are taken from the Jmol repostiory and redistributed under the
GNU Lesser General Public License (see LICENSE.TXT).  The primary goal is simply
to have the g3d classes work.  Accordingly, all classes from org.jmol.g3d are
included along with their dependencies.  Differences from the jmol repository
are primarily related to limiting dependencies to the rest of Jmol or other
libraries.

* Graphics3D is the only implementation of JmolRednererInterface.
* All references and methods related to JmolViewer have been removed.
* All methods, classes and references related to JmolMouse, JmolFile*,
   JmolPopup* and JMolFrame have been removed.
* Methods related to JPEG have been removed.
* All references and methods related to JSObject have been removed.
* All references and methods related to Atom and BondSet have been removed.
* References to J2SRequireImport have been removed.

Many classes and methods that are included are unnecessary, but are included to
minimize the number of changes that need to be made.

The original Jmol source is available for download:
http://jmol.sourceforge.net/download/
https://jmol.svn.sourceforge.net/svnroot/jmol/trunk/Jmol

Also included are classes from Kenji Hiranabe's unofficial Java3d vecmath
package. 
  Copyright (C) 1997,1998,1999 by Kenji Hiranabe.
  http://objectclub.jp/download/vecmath_e