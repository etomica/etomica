package org.jmol.api;


import org.jmol.util.JmolFont;
import org.jmol.util.Point3f;

public interface ApiPlatform {

  /////// Display

  boolean isHeadless();
  
  void convertPointFromScreen(Object display, Point3f ptTemp);

  void getFullScreenDimensions(Object display, int[] widthHeight);
  
  boolean hasFocus(Object display);

  String prompt(String label, String data, String[] list, boolean asButtons);

  void repaint(Object display);

  void requestFocusInWindow(Object display);

  void setCursor(int i, Object display);

  void setTransparentCursor(Object display);

  ///// Font
  
  int fontStringWidth(JmolFont font, Object fontMetrics, String text);

  int getFontAscent(Object fontMetrics);

  int getFontDescent(Object fontMetrics);

  Object getFontMetrics(JmolFont font, Object graphics);

  Object newFont(String fontFace, boolean isBold, boolean isItalic, float fontSize);

  ///// core Image handling
  
  Object allocateRgbImage(int windowWidth, int windowHeight, int[] pBuffer,
                          int windowSize, boolean backgroundTransparent);

  void disposeGraphics(Object graphicForText);

  void drawImage(Object g, Object img, int x, int y, int width, int height);

  int[] drawImageToBuffer(Object gObj, Object imageOffscreen,
                          Object image, int width, int height, int bgcolor);

  void flushImage(Object imagePixelBuffer);

  Object getStaticGraphics(Object image, boolean backgroundTransparent);

  Object getGraphics(Object image);

  int getImageWidth(Object image);

  int getImageHeight(Object image);

  Object newBufferedImage(Object image, int i, int height);

  Object newOffScreenImage(int w, int h);
  
  int[] getTextPixels(String text, JmolFont font3d, Object gObj,
                      Object image, int mapWidth, int height,
                      int ascent);

  ///// Image creation for export (optional for any platform)

  /**
   * can be ignored (return null) if platform cannot save images
   * 
   * @param ret
   * @return     null only if this platform cannot save images
   */
  Object createImage(Object ret);

  /**
   * 
   * @param image
   * @param width
   * @param height
   * @param pixels 
   * @param startRow 
   * @param nRows 
   * @return         pixels
   */
  int[] grabPixels(Object image, int width, int height, 
                   int[] pixels, int startRow, int nRows);

  boolean isSingleThreaded();

  void notifyEndOfRendering();

  public String getDateFormat();

}
