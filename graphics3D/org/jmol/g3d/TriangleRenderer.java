/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2003-2005  Miguel, Jmol Development, www.jmol.org
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

package org.jmol.g3d;


//import org.jmol.util.Logger;
import org.jmol.util.GData;
import org.jmol.util.Point3f;
import org.jmol.util.Point3i;
import org.jmol.util.Rgb16;

/**
 * renders triangles
 *<p>
 * currently only renders flat triangles
 *<p>
 * will probably need performance tuning
 *
 * @author Miguel, miguel@jmol.org
 */
class TriangleRenderer {

  final Graphics3D g3d;
  final LineRenderer line3d;

  private final static int DEFAULT = 64;
  private int[] ax = new int[3], ay = new int[3], az = new int[3];
  private int[] axW = new int[DEFAULT], azW = new int[DEFAULT];
  private int[] axE = new int[DEFAULT], azE = new int[DEFAULT];

  private Rgb16[] rgb16sW, rgb16sE;
  private Rgb16[] rgb16sGouraud;
  
  //private final static boolean VERIFY = false;

  TriangleRenderer(Graphics3D g3d) {
    rgb16sW = new Rgb16[DEFAULT];
    rgb16sE = new Rgb16[DEFAULT];
    for (int i = DEFAULT; --i >= 0;) {
      rgb16sW[i] = new Rgb16();
      rgb16sE[i] = new Rgb16();
    }
    this.g3d = g3d;
    this.line3d = g3d.line3d;
    rgb16sGouraud = new Rgb16[3];
    for (int i = 3; --i >= 0;)
      rgb16sGouraud[i] = new Rgb16();
  }

  Rgb16[] reallocRgb16s(Rgb16[] rgb16s, int n) {
    Rgb16[] t = new Rgb16[n];
    System.arraycopy(rgb16s, 0, t, 0, rgb16s.length);
    for (int i = rgb16s.length; i < n; ++i)
      t[i] = new Rgb16();
    return t;
  }

  final Rgb16 rgb16t1 = new Rgb16();
  final Rgb16 rgb16t2 = new Rgb16();

  void setGouraud(int rgbA, int rgbB, int rgbC) {
    rgb16sGouraud[0].setInt(rgbA);
    rgb16sGouraud[1].setInt(rgbB);
    rgb16sGouraud[2].setInt(rgbC);
  }

  /*===============================================================
   * 2004 05 12 - mth
   * I have been working hard to get the triangles to render
   * correctly when lines are drawn only once.
   * the rules were :
   * a pixel gets drawn when
   * 1. it is to the left of a line
   * 2. it is under a horizontal line
   *
   * this generally worked OK, but failed on small skinny triangles
   * careful reading of Michael Abrash's book
   * Graphics Programming Black Book
   * Chapter 38, The Polygon Primeval, page 714
   * it says:
   *   Narrow wedges and one-pixel-wide polygons will show up spottily
   * I do not understand why this is the case
   * so, the triangle drawing now paints overlapping edges by one pixel
   *
   * note added 10/20/2012 -- Bob Hanson
   * 
   *   These problems were addressed and fixed in 
   *   revision (to Triangle3D) 7086 and 7148 back in 2007.
   *   http://jmol.svn.sourceforge.net/viewvc/jmol/trunk/Jmol/src/org/jmol/g3d/Triangle3D.java?r1=7084&r2=7086&pathrev=17476
   *
   *==============================================================*/

  void drawfillTriangle(int xA, int yA, int zA, int xB, int yB, int zB, int xC,
                        int yC, int zC, boolean useGouraud) {
    ax[0] = xA;
    ax[1] = xB;
    ax[2] = xC;
    ay[0] = yA;
    ay[1] = yB;
    ay[2] = yC;
    az[0] = zA;
    az[1] = zB;
    az[2] = zC;
    fillTriangleB(useGouraud);
  }

  void fillTriangleXYZ(int xScreenA, int yScreenA, int zScreenA, int xScreenB,
                    int yScreenB, int zScreenB, int xScreenC, int yScreenC,
                    int zScreenC, boolean useGouraud) {
    ax[0] = xScreenA;
    ax[1] = xScreenB;
    ax[2] = xScreenC;
    ay[0] = yScreenA;
    ay[1] = yScreenB;
    ay[2] = yScreenC;
    az[0] = zScreenA;
    az[1] = zScreenB;
    az[2] = zScreenC;
    fillTriangleB(useGouraud);
  }

  void fillTriangleP3i(Point3i screenA, Point3i screenB, Point3i screenC,
                    boolean useGouraud) {
    ax[0] = screenA.x;
    ax[1] = screenB.x;
    ax[2] = screenC.x;
    ay[0] = screenA.y;
    ay[1] = screenB.y;
    ay[2] = screenC.y;
    az[0] = screenA.z;
    az[1] = screenB.z;
    az[2] = screenC.z;
    fillTriangleB(useGouraud);
  }

  void fillTriangleP3f(Point3f screenA, Point3f screenB, Point3f screenC,
                    boolean useGouraud) {
    //if (screenA.y > 260)return;
    ax[0] = Math.round(screenA.x);
    ax[1] = Math.round(screenB.x);
    ax[2] = Math.round(screenC.x);
    ay[0] = Math.round(screenA.y);
    ay[1] = Math.round(screenB.y);
    ay[2] = Math.round(screenC.y);
    az[0] = Math.round(screenA.z);
    az[1] = Math.round(screenB.z);
    az[2] = Math.round(screenC.z);
    fillTriangleB(useGouraud);
  }

  void fillTriangleP3if(Point3i screenA, Point3i screenB, Point3i screenC,
                    float factor, boolean useGouraud) {
    ax[0] = screenA.x;
    ax[1] = screenB.x;
    ax[2] = screenC.x;
    ay[0] = screenA.y;
    ay[1] = screenB.y;
    ay[2] = screenC.y;
    az[0] = screenA.z;
    az[1] = screenB.z;
    az[2] = screenC.z;
    adjustVertex(ax, factor);
    adjustVertex(ay, factor);
    adjustVertex(az, factor);
    fillTriangleB(useGouraud);
  }

  private static void adjustVertex(int[] t, float factor) {
    float av = (t[0] + t[1] + t[2]) / 3f;
    for (int i = 0; i < 3; i++)
      t[i] += factor * (av - t[i]);
  }

  private void fillTriangleB(boolean useGouraud) {
    if (az[0] <= 1 || az[1] <= 1 || az[2] <= 1)
      return;
    int cc0 = g3d.clipCode3(ax[0], ay[0], az[0]);
    int cc1 = g3d.clipCode3(ax[1], ay[1], az[1]);
    int cc2 = g3d.clipCode3(ax[2], ay[2], az[2]);
    boolean isClipped = (cc0 | cc1 | cc2) != 0;
    if (isClipped) {
      if ((cc0 & cc1 & cc2) != 0) {
        // all three corners are being clipped on the same dimension
        return;
      }
    }
    int iMinY = 0;
    if (ay[1] < ay[iMinY])
      iMinY = 1;
    if (ay[2] < ay[iMinY])
      iMinY = 2;
    int iMidY = (iMinY + 1) % 3;
    int iMaxY = (iMinY + 2) % 3;
    if (ay[iMidY] > ay[iMaxY]) {
      int t = iMidY;
      iMidY = iMaxY;
      iMaxY = t;
    }
    int yMin = ay[iMinY];
    int yMid = ay[iMidY];
    int yMax = ay[iMaxY];
    int nLines = yMax - yMin + 1;
    if (nLines > g3d.height * 3)
      return;
    if (nLines > axW.length)
      reallocRasterArrays(nLines);
    Rgb16[] gouraudW, gouraudE;
    if (useGouraud) {
      gouraudW = rgb16sW;
      gouraudE = rgb16sE;
    } else {
      gouraudW = gouraudE = null;
    }
    int dyMidMin = yMid - yMin;
    if (dyMidMin == 0) {
      // flat top
      if (ax[iMidY] < ax[iMinY]) {
        int t = iMidY;
        iMidY = iMinY;
        iMinY = t;
      }
      /*  
         min ------- mid
           \         /
          A \       / B
           \ \     / /
              \   /
               max
      */
      generateRaster(nLines, iMinY, iMaxY, axW, azW, 0, gouraudW);
      generateRaster(nLines, iMidY, iMaxY, axE, azE, 0, gouraudE);
    } else if (yMid == yMax) {
      // flat bottom
      if (ax[iMaxY] < ax[iMidY]) {
        int t = iMidY;
        iMidY = iMaxY;
        iMaxY = t;
      }
      /*
       *       min
       *      /   \
       *   A /     \ B
       *  / /       \ \
       *   /         \
       *  mid ------ max
       */
      generateRaster(nLines, iMinY, iMidY, axW, azW, 0, gouraudW);
      generateRaster(nLines, iMinY, iMaxY, axE, azE, 0, gouraudE);
    } else {
      int dxMaxMin = ax[iMaxY] - ax[iMinY];
      int roundFactor;
      roundFactor = GData.roundInt(nLines / 2);
      if (dxMaxMin < 0)
        roundFactor = -roundFactor;
      int axSplit = ax[iMinY] + (dxMaxMin * dyMidMin + roundFactor) / nLines;
      if (axSplit < ax[iMidY]) {
        /*
         *       min
         *      /   \ B
         *   A /     \ \  
         *  / /      mid 
         *   /     C
         *  max   / 
         */

        // Trick is that we need to overlap so as to generate the IDENTICAL
        // raster on each segment, but then we always throw out the FIRST raster
        generateRaster(nLines, iMinY, iMaxY, axW, azW, 0, gouraudW);
        generateRaster(dyMidMin + 1, iMinY, iMidY, axE, azE, 0, gouraudE);
        generateRaster(nLines - dyMidMin, iMidY, iMaxY, axE, azE, dyMidMin,
            gouraudE);

      } else {

        /*
         *       min
         *      /   \ C
         *   A /     \ \  
         *  / /       \ 
         *   /         \
         *  mid         \ 
         *       B->    max
         */

        generateRaster(dyMidMin + 1, iMinY, iMidY, axW, azW, 0, gouraudW);
        generateRaster(nLines - dyMidMin, iMidY, iMaxY, axW, azW, dyMidMin,
            gouraudW);
        generateRaster(nLines, iMinY, iMaxY, axE, azE, 0, gouraudE);
      }
    }
    g3d.setZMargin(5);
    if (useGouraud)
      fillRasterG(yMin, nLines, isClipped, g3d.isPass2 ? 1 : 0);
    else
      fillRaster(yMin, nLines, isClipped, g3d.isPass2 ? 1 : 0);
    g3d.setZMargin(0);
  }

  private void reallocRasterArrays(int n) {
    n = (n + 31) & ~31;
    axW = new int[n];
    azW = new int[n];
    axE = new int[n];
    azE = new int[n];
    rgb16sW = reallocRgb16s(rgb16sW, n);
    rgb16sE = reallocRgb16s(rgb16sE, n);
  }

  private void generateRaster(int dy, int iN, int iS, int[] axRaster,
                              int[] azRaster, int iRaster, Rgb16[] gouraud) {
    int xN = ax[iN], zN = az[iN];
    int xS = ax[iS], zS = az[iS];
    int dx = xS - xN, dz = zS - zN;
    int xCurrent = xN;
    int xIncrement, width, errorTerm;
    if (dx >= 0) {
      xIncrement = 1;
      width = dx;
      errorTerm = 0;
    } else {
      xIncrement = -1;
      width = -dx;
      errorTerm = 1 - dy;
    }
    int zCurrentScaled = (zN << 10) + (1 << 9);
    int roundingFactor = GData.roundInt(dy / 2);
    if (dz < 0)
      roundingFactor = -roundingFactor;
    int zIncrementScaled = ((dz << 10) + roundingFactor) / dy;

    int xMajorIncrement;
    int xMajorError;  
    if (width <= dy) {
      // high-slope
      xMajorIncrement = 0;
      xMajorError = width;
    } else {
      // low-slope
      xMajorIncrement = GData.roundInt(dx / dy);
      xMajorError = width % dy;
    }
    for (int y = 0, i = iRaster; y < dy; zCurrentScaled += zIncrementScaled, ++i, ++y) {
      axRaster[i] = xCurrent;
      azRaster[i] = zCurrentScaled >> 10;
      xCurrent += xMajorIncrement;
      errorTerm += xMajorError;
      if (errorTerm > 0) {
        xCurrent += xIncrement;
        errorTerm -= dy;
      }
    }
    if (gouraud != null) {
      Rgb16 rgb16Base = rgb16t1;
      rgb16Base.setRgb(rgb16sGouraud[iN]);
      Rgb16 rgb16Increment = rgb16t2;
      rgb16Increment.diffDiv(rgb16sGouraud[iS], rgb16Base, dy);
      for (int i = iRaster, iMax = iRaster + dy; i < iMax; ++i)
        gouraud[i].setAndIncrement(rgb16Base, rgb16Increment);
//      if (VERIFY) {
//        Rgb16 north = rgb16sGouraud[iN];
//        Rgb16 generated = gouraud[iRaster];
//        if (north.getArgb() != generated.getArgb()) {
//          if (Logger.debugging) {
//            Logger.debug("north=" + north + "\ngenerated=" + generated);
//          }
//          throw new NullPointerException();
//        }
//      }
    }
  }

  private void fillRaster(int y, int numLines,
                          boolean isClipped, int correction) {
    int i = 0;
    if (y < 0) {
      numLines += y;
      i -= y;
      y = 0;
    }
    if (y + numLines > g3d.height)
      numLines = g3d.height - y;
    if (isClipped) {
      for (; --numLines >= correction; ++y, ++i) {
        int xW = axW[i];
        int pixelCount = axE[i] - xW + 1 - correction;
        if (pixelCount > 0)
          g3d.plotPixelsClippedRaster(pixelCount, xW, y, azW[i], azE[i], null, null);
      }
    } else {
      int xW;
      for (; --numLines >= correction; ++y, ++i) {
        int pixelCount = axE[i] - (xW = axW[i]) + 1 - correction;
        if (correction == 1 && pixelCount < 0) {
          /*
           * The issue here is that some very long, narrow triangles can be skipped
           * altogether because axE < axW.
           * 
           */

          pixelCount = 1;
          xW--;
        }
        if (pixelCount > 0)
          g3d.plotPixelsUnclippedRaster(pixelCount, xW, y, azW[i], azE[i], null, null);
      }
    }
  }

  private void fillRasterG(int y, int numLines,
                          boolean isClipped, int correction) {
    int i = 0;
    if (y < 0) {
      numLines += y;
      i -= y;
      y = 0;
    }
    if (y + numLines > g3d.height)
      numLines = g3d.height - y;
    if (isClipped) {
      for (; --numLines >= correction; ++y, ++i) {
        int xW = axW[i];
        int pixelCount = axE[i] - xW + 1 - correction;
        if (pixelCount > 0)
          g3d.plotPixelsClippedRaster(pixelCount, xW, y, azW[i], azE[i], rgb16sW[i], rgb16sE[i]);
      }
    } else {
      int xW;
      for (; --numLines >= correction; ++y, ++i) {
        int pixelCount = axE[i] - (xW = axW[i]) + 1 - correction;
        if (correction == 1 && pixelCount < 0) {
          /*
           * The issue here is that some very long, narrow triangles can be skipped
           * altogether because axE < axW.
           * 
           */

          pixelCount = 1;
          xW--;
        }
        if (pixelCount > 0)
          g3d.plotPixelsUnclippedRaster(pixelCount, xW, y, azW[i], azE[i], rgb16sW[i], rgb16sE[i]);
      }
    }
  }

}
