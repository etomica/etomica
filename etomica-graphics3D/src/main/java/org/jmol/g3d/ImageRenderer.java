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


class ImageRenderer {

  /**
   * 
   * @param x
   * @param y
   * @param z
   * @param image
   * @param g3d
   * @param jmolRenderer
   * @param antialias  UNUSED
   * @param argbBackground
   * @param width
   * @param height
   */
  static void plotImage(int x, int y, int z, Object image,
                               Graphics3D g3d,
                               Graphics3D jmolRenderer,
                               boolean antialias, int argbBackground,
                               int width, int height) {
    boolean isBackground = (x == Integer.MIN_VALUE);
    int bgcolor = (isBackground ? g3d.bgcolor : argbBackground);
    /*
     *  this was for transparent background, which we have disabled, I think, in Jmol 12
    boolean haveTranslucent = false;
    PixelGrabber pg1 = new PixelGrabber(image, 0, 0, width0, height0, true);
    if (pg1.getColorModel().hasAlpha())
      try {
        pg1.grabPixels();
        int[] buffer = (int[]) pg1.getPixels();
        for (int i = 0; i < buffer.length; i++)
          if ((buffer[i] & 0xFF00000) != 0xFF000000) {
            haveTranslucent = true;
            break;
          }
        System.out.println(buffer.length + " " + haveTranslucent + " "
            + pg1.getColorModel().hasAlpha());
      } catch (InterruptedException e) {
        // impossible?
        return;
      }
      */
    if (isBackground) {
      x = 0;
      z = Integer.MAX_VALUE - 1;
      width = g3d.width;
      height = g3d.height;
    }
    if (x + width <= 0 || x >= g3d.width || y + height <= 0 || y >= g3d.height)
      return;
    Object g;
    /**
     * @j2sNative
     * 
     * g = null;
     *  
     */
    {
      g = g3d.platform.getGraphicsForTextOrImage(width, height);
    } 
    int[] buffer = g3d.apiPlatform.drawImageToBuffer(g, g3d.platform.offscreenImage, image, width, height, 
        isBackground ? bgcolor : 0);
    if (buffer == null)
      return; // not supported on this platform (yet)
/*    
    int n = 0;
    for (int i = 0; i < buffer.length; i++) {
      if ((buffer[i] & 0xFF000000) != 0xFF000000) {
        //System.out.println("testing " + i + " " + buffer[i]);
        n++;
      }
    }
    System.out.println(n + " transparent argbBackground=" + argbBackground);
*/
    if (jmolRenderer != null
        || (x < 0 || x + width > g3d.width || y < 0 || y + height > g3d.height))
      plotImageClipped(x, y, z, g3d, jmolRenderer, width, height, buffer,
          bgcolor);
    else
      plotImageUnClipped(x, y, z, g3d, width, height, buffer, bgcolor);
  }

  private static void plotImageClipped(int x, int y, int z, Graphics3D g3d,
                                       Graphics3D jmolRenderer,
                                       int width, int height,
                                       int[] buffer, int bgargb) {
    if (jmolRenderer == null)
      jmolRenderer = g3d;
    for (int i = 0, offset = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        int argb = buffer[offset++];
          jmolRenderer.plotImagePixel(argb, x + j, y + i, z, 8, bgargb);
      }
    }
  }

  /**
   * @param x 
   * @param y 
   * @param z 
   * @param g3d 
   * @param textWidth 
   * @param textHeight 
   * @param buffer 
   * @param bgargb  used in transparent backgrounds?  
   */
  private static void plotImageUnClipped(int x, int y, int z, Graphics3D g3d,
                                         int textWidth, int textHeight,
                                         int[] buffer, int bgargb) {
    int[] zbuf = g3d.zbuf;
    int renderWidth = g3d.width;
    int pbufOffset = y * renderWidth + x;
    int i = 0;
    int j = 0;
    int offset = 0;
    while (i < textHeight) {
      while (j < textWidth) {
        if (z < zbuf[pbufOffset]) {
          int argb = buffer[offset];
          //if (argb != bgargb) && (argb & 0xFF000000) == 0xFF000000)
            g3d.addPixel(pbufOffset, z, argb);
          //else if (argb == 0 && bgargb != 0)
            //g3d.addPixel(pbufOffset, z, bgargb);
        }
        ++offset;
        ++j;
        ++pbufOffset;
      }
      ++i;
      j -= textWidth;
      pbufOffset += (renderWidth - textWidth);
    }
  }

}
