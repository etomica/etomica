/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2003-2005  Miguel, Jmol Development, www.jmol.org
 *
 * Contact: miguel@jmol.org
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
package g3dsys.images;

import org.jmol.g3d.Graphics3D;

/**
 *<p>
 * Implements ellipse drawing/filling routines, based on Circl3D
 *</p>
 *
 * @author Miguel, miguel@jmol.org
 * @author Andrew Schultz, ajs42@buffalo.edu
 */
public class Ellipse3D {

  Graphics3D g3d;

  public Ellipse3D(Graphics3D g3d) {
    this.g3d = g3d;
  }

  int xCenter, yCenter, zCenter;
  int sizeCorrection;
  
  void plotCircleCenteredClipped(int xCenter, int yCenter, int zCenter,
                                 int xDiameter, float aspect) {
    // halo only -- simple window clip
    {
      int r = (xDiameter + 1) >> 1;
      if (xCenter + r < 0 || xCenter - r >= g3d.getRenderWidth()||
          yCenter + r < 0 || yCenter - r >= g3d.getRenderHeight())
        return;
    }
    int r = xDiameter / 2;
    this.sizeCorrection = 1 - (xDiameter & 1);
    this.xCenter = xCenter;
    this.yCenter = yCenter;
    this.zCenter = zCenter;
    int x = r;
    int y = 0;
    int xChange = 1 - 2*r;
    float yChange = aspect;
    int intYChange = 1;
    float radiusError = 0;
    int intRadiusError = 0;
    while (true) {
      plot8CircleCenteredClipped(x, y);
      ++y;
      radiusError += yChange;
      intRadiusError += intYChange;
      yChange += 2*aspect;
      intYChange += 2;
      
      while (2*radiusError + xChange > 0) {
        --x;
        radiusError += xChange;
        intRadiusError += xChange;
        xChange += 2;
        if (2*radiusError + xChange > 0) {
            plot8CircleCenteredClipped(x, y);
        }
        if (x == 0) return;
      }
    }
  }

  private void plot8CircleCenteredClipped(int dx, int dy) {
    g3d.drawPixel(xCenter+dx-sizeCorrection,
                          yCenter+dy-sizeCorrection, zCenter);
    g3d.drawPixel(xCenter+dx-sizeCorrection, yCenter-dy, zCenter);
    g3d.drawPixel(xCenter-dx, yCenter+dy-sizeCorrection, zCenter);
    g3d.drawPixel(xCenter-dx, yCenter-dy, zCenter);
  }
}

