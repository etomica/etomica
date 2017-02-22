/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2006  Miguel, Jmol Development, www.jmol.org
 *
 * Contact: jmol-developers@lists.sf.net
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package org.jmol.util;

/**
 * the Point3fi class allows storage of critical information involving
 * an atom, picked bond, or measurement point, including: 
 * 
 * xyz position 
 * screen position
 * screen radius (-1 for a simple point)
 * index (for atoms or for an associated bond that has be picked)
 * associated modelIndex (for measurement points)
 * 
 */
public class Point3fi extends Point3f {
  public int index;
  public int screenX;
  public int screenY;
  public int screenZ;
  public short screenDiameter = -1;
  public short modelIndex = -1;

  public static void set2(Point3f p3f, Point3i p3i) {
    p3f.x = p3i.x;
    p3f.y = p3i.y;
    p3f.z = p3i.z;
  }

}
