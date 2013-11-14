package org.jmol.util;


import java.awt.Image;


import org.jmol.api.ApiPlatform;
import org.jmol.constant.EnumStereoMode;
import org.jmol.g3d.Graphics3D;

public class GData {

  public ApiPlatform apiPlatform;

  protected int windowWidth, windowHeight;
  protected int displayMinX, displayMaxX, displayMinY, displayMaxY;
  protected boolean antialiasThisFrame;
  protected boolean antialiasEnabled;

  protected boolean inGreyscaleMode;

  protected short[] changeableColixMap = new short[16];
  
  protected Object backgroundImage;
  
  protected int newWindowWidth, newWindowHeight;
  protected boolean newAntialiasing;

  public int bgcolor;  
  public int xLast, yLast;
  public int slab, depth;
  public int width, height;
  public int zSlab, zDepth;
  public int zShadePower = 3;

  protected short colixCurrent;
  protected int argbCurrent;

  public int bufferSize;

  public final static byte ENDCAPS_NONE = 0;
  public final static byte ENDCAPS_OPEN = 1;
  public final static byte ENDCAPS_FLAT = 2;
  public final static byte ENDCAPS_SPHERICAL = 3;
  public final static byte ENDCAPS_OPENEND = 4;
  
  public GData() {
    
  }
  
  public void initialize(ApiPlatform apiPlatform) {
    this.apiPlatform = apiPlatform;
  }

  /**
   * clipping from the front and the back
   *<p>
   * the plane is defined as a percentage from the back of the image to the
   * front
   *<p>
   * for depth values:
   * <ul>
   * <li>0 means 100% is shown
   * <li>25 means the back 25% is <i>not</i> shown
   * <li>50 means the back half is <i>not</i> shown
   * <li>100 means that nothing is shown
   * </ul>
   *<p>
   * 
   * @param depthValue
   *        rear clipping percentage [0,100]
   */
  public void setDepth(int depthValue) {
    depth = depthValue < 0 ? 0 : depthValue;
  }

  /**
   * clipping from the front and the back
   *<p>
   * the plane is defined as a percentage from the back of the image to the
   * front
   *<p>
   * For slab values:
   * <ul>
   * <li>100 means 100% is shown
   * <li>75 means the back 75% is shown
   * <li>50 means the back half is shown
   * <li>0 means that nothing is shown
   * </ul>
   *<p>
   * 
   * @param slabValue
   *        front clipping percentage [0,100]
   */
  public void setSlab(int slabValue) {
    slab = slabValue < 0 ? 0 : slabValue;
  }

  public int zShadeR, zShadeG, zShadeB;

  protected Object graphicsForMetrics;

  public final static int EXPORT_RAYTRACER = 2;

  public final static int EXPORT_CARTESIAN = 1;

  public final static int EXPORT_NOT = 0;

  /**
   * @param zShade
   *        whether to shade along z front to back
   * @param zSlab
   *        for zShade
   * @param zDepth
   *        for zShade
   * @param zPower 
   */
  public void setZShade(boolean zShade, int zSlab, int zDepth, int zPower) {
    if (zShade) {
      zShadeR = bgcolor & 0xFF;
      zShadeG = (bgcolor & 0xFF00) >> 8;
      zShadeB = (bgcolor & 0xFF0000) >> 16;
      this.zSlab = zSlab < 0 ? 0 : zSlab;
      this.zDepth = zDepth < 0 ? 0 : zDepth;
      this.zShadePower = zPower;
    }
  }

  /**
   * gets g3d width
   * 
   * @return width pixel count;
   */
  public int getRenderWidth() {
    return width;
  }

  /**
   * gets g3d height
   * 
   * @return height pixel count
   */
  public int getRenderHeight() {
    return height;
  }

  /**
   * gets g3d slab
   * 
   * @return slab
   */
  public int getSlab() {
    return slab;
  }

  /**
   * gets g3d depth
   * 
   * @return depth
   */
  public int getDepth() {
    return depth;
  }


  /**
   * is full scene / oversampling antialiasing GENERALLY in effect
   *
   * @return the answer
   */
  public boolean isDisplayAntialiased() {
    return antialiasEnabled;
  }

  /**
   * is full scene / oversampling antialiasing in effect
   * 
   * @return the answer
   */
  public boolean isAntialiased() {
    return antialiasThisFrame;
  }

  public short getChangeableColix(short id, int argb) {
    if (id >= changeableColixMap.length)
      changeableColixMap = ArrayUtil.arrayCopyShort(changeableColixMap, id + 16);
    if (changeableColixMap[id] == 0)
      changeableColixMap[id] = Colix.getColix(argb);
    return (short) (id | Colix.CHANGEABLE_MASK);
  }

  public void changeColixArgb(short id, int argb) {
    if (id < changeableColixMap.length && changeableColixMap[id] != 0)
      changeableColixMap[id] = Colix.getColix(argb);
  }

  public short[] getBgColixes(short[] bgcolixes) {
    return bgcolixes;
  }

  public int getColorArgbOrGray(short colix) {
    if (colix < 0)
      colix = changeableColixMap[colix & Colix.UNMASK_CHANGEABLE_TRANSLUCENT];
    return (inGreyscaleMode ? Colix.getArgbGreyscale(colix) : Colix.getArgb(colix));
  }

  public int[] getShades(short colix) {
    if (colix < 0)
      colix = changeableColixMap[colix & Colix.UNMASK_CHANGEABLE_TRANSLUCENT];
    return (inGreyscaleMode ? Colix.getShadesGreyscale(colix) : Colix
        .getShades(colix));
  }

  /**
   * controls greyscale rendering
   * 
   * @param greyscaleMode
   *        Flag for greyscale rendering
   */
  public void setGreyscaleMode(boolean greyscaleMode) {
    this.inGreyscaleMode = greyscaleMode;
  }

  public static int getSpecularPower() {
    return Shader.specularPower;
  }

  /**
   * fractional distance to white for specular dot
   * 
   * @param val
   */
  public synchronized static void setSpecularPower(int val) {
    if (val < 0) {
      setSpecularExponent(-val);
      return;
    }
    if (Shader.specularPower == val)
      return;
    Shader.specularPower = val;
    Shader.intenseFraction = val / 100f;
    flushCaches();
  }

  public static int getSpecularPercent() {
    return Shader.specularPercent;
  }

  /**
   * sf in I = df * (N dot L) + sf * (R dot V)^p not a percent of anything,
   * really
   * 
   * @param val
   */
  public synchronized static void setSpecularPercent(int val) {
    if (Shader.specularPercent == val)
      return;
    Shader.specularPercent = val;
    Shader.specularFactor = val / 100f;
    flushCaches();
  }

  public static int getSpecularExponent() {
    return Shader.specularExponent;
  }

  /**
   * log_2(p) in I = df * (N dot L) + sf * (R dot V)^p for faster calculation of
   * shades
   * 
   * @param val
   */
  public synchronized static void setSpecularExponent(int val) {
    if (Shader.specularExponent == val)
      return;
    Shader.specularExponent = val;
    Shader.phongExponent = (int) Math.pow(2, val);
    Shader.usePhongExponent = false;
    flushCaches();
  }

  public static int getPhongExponent() {
    return Shader.phongExponent;
  }

  /**
   * p in I = df * (N dot L) + sf * (R dot V)^p
   * 
   * @param val
   */
  public synchronized static void setPhongExponent(int val) {
    if (Shader.phongExponent == val && Shader.usePhongExponent)
      return;
    Shader.phongExponent = val;
    float x = (float) (Math.log(val) / Math.log(2));
    Shader.usePhongExponent = (x != (int) x);
    if (!Shader.usePhongExponent)
      Shader.specularExponent = (int) x;
    flushCaches();
  }

  public static int getDiffusePercent() {
    return Shader.diffusePercent;
  }

  /**
   * df in I = df * (N dot L) + sf * (R dot V)^p
   * 
   * @param val
   */
  public synchronized static void setDiffusePercent(int val) {
    if (Shader.diffusePercent == val)
      return;
    Shader.diffusePercent = val;
    Shader.diffuseFactor = val / 100f;
    flushCaches();
  }

  public static int getAmbientPercent() {
    return Shader.ambientPercent;
  }

  /**
   * fractional distance from black for ambient color
   * 
   * @param val
   */
  public synchronized static void setAmbientPercent(int val) {
    if (Shader.ambientPercent == val)
      return;
    Shader.ambientPercent = val;
    Shader.ambientFraction = val / 100f;
    flushCaches();
  }

  public static boolean getSpecular() {
    return Shader.specularOn;
  }

  public synchronized static void setSpecular(boolean val) {
    if (Shader.specularOn == val)
      return;
    Shader.specularOn = val;
    flushCaches();
  }

  private static void flushCaches() {
    Colix.flushShades();
    Shader.flushSphereCache();
  }

  public static Point3f getLightSource() {
    return Point3f.new3(Shader.xLight, Shader.yLight, Shader.zLight);
  }

  public boolean isClipped3(int x, int y, int z) {
    // this is the one that could be augmented with slabPlane
    return (x < 0 || x >= width || y < 0 || y >= height || z < slab || z > depth);
  }

  public boolean isClipped(int x, int y) {
    return (x < 0 || x >= width || y < 0 || y >= height);
  }

  public boolean isInDisplayRange(int x, int y) {
    return (x >= displayMinX && x < displayMaxX && y >= displayMinY && y < displayMaxY);
  }

  public boolean isClippedXY(int diameter, int x, int y) {
    int r = (diameter + 1) >> 1;
    return (x < -r || x >= width + r || y < -r || y >= height + r);
  }

  public boolean isClippedZ(int z) {
    return (z != Integer.MIN_VALUE && (z < slab || z > depth));
  }

  final public static int yGT = 1;
  final public static int yLT = 2;
  final public static int xGT = 4;
  final public static int xLT = 8;
  final public static int zGT = 16;
  final public static int zLT = 32;

  public int clipCode3(int x, int y, int z) {
    int code = 0;
    if (x < 0)
      code |= xLT;
    else if (x >= width)
      code |= xGT;
    if (y < 0)
      code |= yLT;
    else if (y >= height)
      code |= yGT;
    if (z < slab)
      code |= zLT;
    else if (z > depth) // note that this is .GT., not .GE.
      code |= zGT;

    return code;
  }

  public int clipCode(int z) {
    int code = 0;
    if (z < slab)
      code |= zLT;
    else if (z > depth) // note that this is .GT., not .GE.
      code |= zGT;
    return code;
  }

  /* ***************************************************************
   * fontID stuff
   * a fontID is a byte that contains the size + the face + the style
   * ***************************************************************/

  public JmolFont getFont3D(float fontSize) {
    return JmolFont.createFont3D(JmolFont.FONT_FACE_SANS, JmolFont.FONT_STYLE_PLAIN,
        fontSize, fontSize, apiPlatform, graphicsForMetrics);
  }

  public JmolFont getFont3DFS(String fontFace, float fontSize) {
    return JmolFont.createFont3D(JmolFont.getFontFaceID(fontFace),
        JmolFont.FONT_STYLE_PLAIN, fontSize, fontSize, apiPlatform, graphicsForMetrics);
  }

  public byte getFontFidFS(String fontFace, float fontSize) {
    return getFont3DFSS(fontFace, "Bold", fontSize).fid;
  }

  public JmolFont getFont3DFSS(String fontFace, String fontStyle, float fontSize) {
    int iStyle = JmolFont.getFontStyleID(fontStyle);
    if (iStyle < 0)
      iStyle = 0;
    return JmolFont.createFont3D(JmolFont.getFontFaceID(fontFace), iStyle, fontSize,
        fontSize, apiPlatform, graphicsForMetrics);
  }

  public JmolFont getFont3DScaled(JmolFont font, float scale) {
    // TODO: problem here is that we are assigning a bold font, then not DEassigning it
    float newScale = font.fontSizeNominal * scale;
    return (newScale == font.fontSize ? font : JmolFont.createFont3D(
        font.idFontFace, font.idFontStyle, newScale, font.fontSizeNominal, apiPlatform, graphicsForMetrics));
  }

  public byte getFontFid(float fontSize) {
    return getFont3D(fontSize).fid;
  }

  // {"Plain", "Bold", "Italic", "BoldItalic"};
  public static int getFontStyleID(String fontStyle) {
    return JmolFont.getFontStyleID(fontStyle);
  }

  /**
   * @param TF  
   */
  public void setBackgroundTransparent(boolean TF) {
  }

  /**
   * sets background color to the specified argb value
   *
   * @param argb an argb value with alpha channel
   */
  public void setBackgroundArgb(int argb) {
    bgcolor = argb;
    // background of Jmol transparent in front of certain applications (VLC Player)
    // when background [0,0,1]. 
  }

  public void setBackgroundImage(Object image) {
    backgroundImage = image;
  }

  public void setWindowParameters(int width, int height, boolean antialias) {
    newWindowWidth = width;
    newWindowHeight = height;
    newAntialiasing = antialias;    
  }

  public void setNewWindowParametersForExport() {
    windowWidth = newWindowWidth;
    windowHeight = newWindowHeight;
    setWidthHeight(false);
  }

  protected void setWidthHeight(boolean isAntialiased) {
    width = windowWidth;
    height = windowHeight;
    if (isAntialiased) {
      width <<= 1;
      height <<= 1;
    }
    xLast = width - 1;
    yLast = height - 1;
    displayMinX = -(width >> 1);
    displayMaxX = width - displayMinX;
    displayMinY = -(height >> 1);
    displayMaxY = height - displayMinY;
    bufferSize = width * height;
  }

  /**
   * @param stereoRotationMatrix  
   */
  public void beginRendering(Matrix3f stereoRotationMatrix) {
  }

  public void endRendering() {
  }

  public void snapshotAnaglyphChannelBytes() {
  }

  public Object getScreenImage() {
    return null;
  }

  public void releaseScreenImage() {
  }

  /**
   * @param stereoMode  
   * @param stereoColors 
   */
  public void applyAnaglygh(EnumStereoMode stereoMode, int[] stereoColors) {
  }

  /**
   * @param antialias  
   * @return true if need a second (translucent) pass
   */
  public boolean setPass2(boolean antialias) {
    return false;
  }

  public void destroy() {
  }

  public void clearFontCache() {
  }

  /**
   * @param x  
   * @param y 
   * @param z 
   * @param image 
   * @param jmolRenderer
   * @param bgcolix
   * @param width
   * @param height
   *  
   */
  public void plotImage(int x, int y, int z, Image image,
                        Graphics3D jmolRenderer, short bgcolix,
                        int width, int height) {
  }

  /**
   * @param x  
   * @param y 
   * @param z 
   * @param colorArgbOrGray
   * @param bgColor TODO
   * @param text
   * @param font3d 
   * @param jmolRenderer
   *  
   */
  public void plotText(int x, int y, int z, int colorArgbOrGray, int bgColor,
                       String text, JmolFont font3d, Graphics3D jmolRenderer) {
  }

  /**
   * @param jmolRenderer  
   */
  public void renderBackground(Graphics3D jmolRenderer) {
  }

  public JmolFont getFont3DCurrent() {
    return null;
  }

  /**
   * @param font3d  
   */
  public void setFont(JmolFont font3d) {
  }

  /**
   * @param fid  
   */
  public void setFontFid(byte fid) {
  }

  /**
   * @param color  
   */
  public void setColor(int color) {
    argbCurrent = color;
  }

  public boolean isPass2;
  
  public boolean isPass2() {
    return isPass2;
  }  

  /**
   * @param colix  
   * @return TRUE if correct pass (translucent or opaque)
   */
  public boolean setColix(short colix) {    
    return true;
  }

  /**
   * @param normix  
   * @return true if front
   */
  public boolean isDirectedTowardsCamera(short normix) {
    return true;
  }

  public Vector3f[] getTransformedVertexVectors() {
    return null;
  }

  /**
   * @param pointA  
   * @param pointB  
   * @param pointC  
   */
  public void setNoisySurfaceShade(Point3i pointA, Point3i pointB,
                                   Point3i pointC) {
  }

  /**
   * JavaScript won't really have an integer here after integer division.
   * So we need to round it to the integer between it and zero. 
   * 
   * @param a
   * @return number closest to zero
   */
  public static int roundInt(int a) {
    /**
     * @z2sNative
     *
     *  return a < 0 ? Math.ceil(a) : Math.floor(a);
     *  
     */
    {
      return a;
    }
  }

  public void clear() {
    // only in Graphics3D
  }

  public void renderAllStrings(Object jmolRenderer) {
    // only in Graphics3D
  }
}
