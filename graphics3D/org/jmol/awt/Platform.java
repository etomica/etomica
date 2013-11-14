package org.jmol.awt;

import java.awt.Container;
import java.awt.Frame;
import java.awt.GraphicsEnvironment;
import java.awt.Window;
import java.text.SimpleDateFormat;
import java.util.Date;

import javax.swing.JDialog;

import org.jmol.api.ApiPlatform;
import org.jmol.util.JmolFont;
import org.jmol.util.Point3f;

public class Platform implements ApiPlatform {

  ///// Display 

  public void convertPointFromScreen(Object display, Point3f ptTemp) {
    Display.convertPointFromScreen(display, ptTemp);
  }

  public void getFullScreenDimensions(Object display, int[] widthHeight) {
    Display.getFullScreenDimensions(display, widthHeight);        
  }

  public boolean hasFocus(Object display) {
    return Display.hasFocus(display);
  }

  public String prompt(String label, String data, String[] list,
                       boolean asButtons) {
    return Display.prompt(label, data, list, asButtons);
  }

  public void requestFocusInWindow(Object display) {
    Display.requestFocusInWindow(display);
  }

  public void repaint(Object display) {
    Display.repaint(display);
  }

  public void setTransparentCursor(Object display) {
    Display.setTransparentCursor(display);
  }

  public void setCursor(int c, Object display) {
    Display.setCursor(c, display);
  }

  ////// Image 

  public Object allocateRgbImage(int windowWidth, int windowHeight,
                                 int[] pBuffer, int windowSize,
                                 boolean backgroundTransparent) {
    return Image.allocateRgbImage(windowWidth, windowHeight, pBuffer, windowSize, backgroundTransparent);
  }

  /**
   * could be byte[] (from ZIP file) or String (local file name) or URL
   * @param data 
   * @return image object
   * 
   */
  public Object createImage(Object data) {
    return Image.createImage(data);
  }

  public void disposeGraphics(Object gOffscreen) {
    Image.disposeGraphics(gOffscreen);
  }

  public void drawImage(Object g, Object img, int x, int y, int width, int height) {
    Image.drawImage(g, img, x, y, width, height);
  }

  public int[] grabPixels(Object imageobj, int width, int height, int[] pixels, int startRow, int nRows) {
    return Image.grabPixels(imageobj, width, height, pixels, startRow, nRows); 
  }

  public int[] drawImageToBuffer(Object gOffscreen, Object imageOffscreen,
                                 Object imageobj, int width, int height, int bgcolor) {
    return Image.drawImageToBuffer(gOffscreen, imageOffscreen, imageobj, width, height, bgcolor);
  }

  public int[] getTextPixels(String text, JmolFont font3d, Object gObj,
                             Object image, int width, int height, int ascent) {
    return Image.getTextPixels(text, font3d, gObj, image, width, height, ascent);
  }

  public void flushImage(Object imagePixelBuffer) {
    Image.flush(imagePixelBuffer);
  }

  public Object getGraphics(Object image) {
    return Image.getGraphics(image);
  }

  public int getImageHeight(Object image) {
    return Image.getHeight(image);
  }

  public int getImageWidth(Object image) {
    return Image.getWidth(image);
  }

  public Object getStaticGraphics(Object image, boolean backgroundTransparent) {
    return Image.getStaticGraphics(image, backgroundTransparent);
  }

  public Object newBufferedImage(Object image, int w, int h) {
    return Image.newBufferedImage(image, w, h);
  }

  public Object newOffScreenImage(int w, int h) {
    return Image.newBufferedImage(w, h);
  }

  ///// FONT
  
  public int fontStringWidth(JmolFont font, Object fontMetrics, String text) {
    return Font.stringWidth(fontMetrics, text);
  }

  public int getFontAscent(Object fontMetrics) {
    return Font.getAscent(fontMetrics);
  }

  public int getFontDescent(Object fontMetrics) {
    return Font.getDescent(fontMetrics);
  }

  public Object getFontMetrics(JmolFont font, Object graphics) {
    return Font.getFontMetrics(font, graphics);
  }

  public Object newFont(String fontFace, boolean isBold, boolean isItalic, float fontSize) {
    return Font.newFont(fontFace, isBold, isItalic, fontSize);
  }

  /// misc

  public boolean isHeadless() {
    return GraphicsEnvironment.isHeadless();
  }

  public boolean isSingleThreaded() {
    return false;
  }

  public void notifyEndOfRendering() {
    // N/A
  }

  /**
   * @param p 
   * @return The hosting frame or JDialog.
   */
  static public Window getWindow(Container p) {
    while (p != null) {
      if (p instanceof Frame)
        return (Frame) p;
      else if (p instanceof JDialog)
        return (JDialog) p;
      p = p.getParent();
    }
    return null;
  }

  public String getDateFormat() {
    return (new SimpleDateFormat("EEE, d MMM yyyy HH:mm:ss Z"))
        .format(new Date());
  }


}
