/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2011  The Jmol Development Team
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
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301, USA.
 */

package org.jmol.g3d;

class PixelatorShaded extends Pixelator {

  /**
   * @param graphics3d
   */
  PixelatorShaded(Graphics3D graphics3d) {
    super(graphics3d);
  }

  @Override
  void addPixel(int offset, int z, int p) {
    if (z > g.zDepth)
      return;
    if (z <= g.zDepth && z >= g.zSlab) {
      int pR = p & 0xFF;
      int pG = (p & 0xFF00) >> 8;
      int pB = (p & 0xFF0000) >> 16;
      int pA = (p & 0xFF000000);
      float f = (float)(g.zDepth - z) / (g.zDepth - g.zSlab);
      if (g.zShadePower > 1) {
        for (int i = 0; i < g.zShadePower; i++)
          f *= f;
      }
      pR = g.zShadeR + (int) (f * (pR - g.zShadeR));
      pG = g.zShadeG + (int) (f * (pG - g.zShadeG));
      pB = g.zShadeB + (int) (f * (pB - g.zShadeB));        
      p = (pB << 16) | (pG << 8) | pR | pA;
    }
    super.addPixel(offset, z, p);
  }
}