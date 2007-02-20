/* Class gutted for g3dsys.
 * 
 * $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2003-2005  The Jmol Development Team
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
package g3dsys.control;

import javax.vecmath.AxisAngle4f;
import javax.vecmath.Matrix3f;
import javax.vecmath.Matrix4f;
import javax.vecmath.Point3f;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3f;

class TransformManager {

  G3DSys g3dsys;
  
  TransformManager(G3DSys gsys) {
    g3dsys = gsys;
  }

  /* ***************************************************************
   * GENERAL METHODS
   ***************************************************************/

  /**
   * Return view to home position; remove rotation, translation, and zoom
   */
  void homePosition() {
    // reset
    setDefaultRotation();
    setRotationCenterAndRadiusXYZ(null, true);
    translateCenterTo(0, 0);
    matrixRotate.setIdentity(); // no rotations
    setZoomEnabled(true);
    zoomToPercent(100);
    scaleFitToScreen();
    //navigationCenter.set(fixedRotationCenter);
  }

  void clear() {
    fixedRotationCenter.set(0, 0, 0);
  }
  
  final static float twoPI = (float) (2 * Math.PI);
  
  final Point3f fixedRotationCenter = new Point3f(0, 0, 0);
  float rotationRadius;
  Point3f rotationCenterDefault;
  float rotationRadiusDefault;
  
  /* ***************************************************************
   * ROTATIONS
   ***************************************************************/

  // this matrix only holds rotations ... no translations
  // however, it cannot be a Matrix3f because we need to multiply it by
  // a matrix4f which contains translations
  //    ...huh?
  private final Matrix3f matrixRotate = new Matrix3f();
  private final Matrix3f matrixTemp3 = new Matrix3f();

  final AxisAngle4f axisangleT = new AxisAngle4f();
  final Point3f pointT = new Point3f();
  final Point3f pointT2 = new Point3f();
  final static float radiansPerDegree = (float) (2 * Math.PI / 360);
  final static float degreesPerRadian = (float) (360 / (2 * Math.PI));

  final static int MAXIMUM_ZOOM_PERCENTAGE = 200000;
  final static int MAXIMUM_ZOOM_PERSPECTIVE_DEPTH = 10000;

  private void setFixedRotationCenter(Point3f center) {
    if (center == null)
      return;
    fixedRotationCenter.set(center);
  }

  void setRotationPointXY(Point3f center) {
    Point3i newCenterScreen = transformPoint(center);
    translateCenterTo(newCenterScreen.x, newCenterScreen.y);
  }

  Vector3f rotationAxis = new Vector3f();
  float rotationRate = 0;

  float setRotateFixed(Point3f center, Vector3f axis, float degrees) {
    setFixedRotationCenter(center);
    rotationAxis.set(axis);
    float radians = degrees * radiansPerDegree;
    rotationRate = degrees;
    return radians;
  }

  void rotateXYBy(int xDelta, int yDelta) {
    // from mouse action
    float f = 1;
    rotateXRadians(xDelta * radiansPerDegree * f);
    rotateYRadians(yDelta * radiansPerDegree * f);
  }

  void rotateZBy(int zDelta) {
    rotateZRadians((float) Math.PI * zDelta / 180);
  }

  void rotateFront() {
    matrixRotate.setIdentity();
  }

  void rotateToX(float angleRadians) {
    matrixRotate.rotX(angleRadians);
  }

  void rotateToY(float angleRadians) {
    matrixRotate.rotY(angleRadians);
  }

  void rotateToZ(float angleRadians) {
    matrixRotate.rotZ(angleRadians);
  }

  synchronized void rotateXRadians(float angleRadians) {
    matrixTemp3.rotX(angleRadians);
    matrixRotate.mul(matrixTemp3, matrixRotate);
  }

  synchronized void rotateYRadians(float angleRadians) {
    if (axesOrientationRasmol)
      angleRadians = -angleRadians;
    matrixTemp3.rotY(angleRadians);
    matrixRotate.mul(matrixTemp3, matrixRotate);
  }

  synchronized void rotateZRadians(float angleRadians) {
    if (axesOrientationRasmol)
      angleRadians = -angleRadians;
    matrixTemp3.rotZ(angleRadians);
    matrixRotate.mul(matrixTemp3, matrixRotate);
  }

  void rotateTo(float x, float y, float z, float degrees) {
    //unused
    if (degrees < .01 && degrees > -.01) {
      matrixRotate.setIdentity();
    } else {
      axisangleT.set(x, y, z, degrees * radiansPerDegree);
      matrixRotate.set(axisangleT);
    }
  }

  /* ***************************************************************
   * TRANSLATIONS
   ****************************************************************/
  float xFixedTranslation;
  float yFixedTranslation;

  void translateXYBy(int xDelta, int yDelta) {
    // mouse action only
    xFixedTranslation += xDelta;
    yFixedTranslation += yDelta;
  }

  void translateToXPercent(float percent) {
    xFixedTranslation = (width / 2) + width * percent / 100;
  }

  void translateToYPercent(float percent) {
    yFixedTranslation = (height / 2) + height * percent / 100;
  }

  void translateToZPercent(float percent) {
    // FIXME who knows what this should be? some type of zoom?
  }

  float getTranslationXPercent() {
    return (xFixedTranslation - width / 2) * 100 / width;
  }

  float getTranslationYPercent() {
    return (yFixedTranslation - height / 2) * 100 / height;
  }

  float getTranslationZPercent() {
    return 0;
  }

  void translateCenterTo(int x, int y) {
    xFixedTranslation = x;
    yFixedTranslation = y;
  }

  Matrix3f getMatrixRotate() {
    return matrixRotate;
  }

  void setRotation(Matrix3f matrixRotation) {
    matrixRotate.set(matrixRotation);
  }

  void getRotation(Matrix3f matrixRotation) {
    // hmm ... I suppose that there could be a race condiditon here
    // if matrixRotate is being modified while this is called
    matrixRotation.set(matrixRotate);
  }

  /* ***************************************************************
   * ZOOM
   ****************************************************************/
  boolean zoomEnabled = true;
  // zoomPercent is the current displayed zoom value
  float zoomPercent = 100;
  // zoomPercentSetting is the current setting of zoom
  // if zoom is not enabled then the two values will be different
  float zoomPercentSetting = 100;

  void zoomBy(int pixels) {
    if (pixels > 20)
      pixels = 20;
    else if (pixels < -20)
      pixels = -20;
    float deltaPercent = pixels * zoomPercentSetting / 50;
    if (deltaPercent == 0)
      deltaPercent = (pixels > 0 ? 1 : (deltaPercent < 0 ? -1 : 0));
    float percent = deltaPercent + zoomPercentSetting;
    zoomToPercent(percent);
  }

  int getZoomPercent() {
    return (int) zoomPercent;
  }

  float getZoomPercentFloat() {
    return zoomPercent;
  }

  float getZoomPercentSetting() {
    return zoomPercentSetting;
  }

  void zoomToPercent(float percentZoom) {
    zoomPercentSetting = percentZoom;
    calcScale("zoomToPercent");
  }

  void zoomByPercent(float percentZoom) {
    float delta = percentZoom * zoomPercentSetting / 100;
    if (delta == 0)
      delta = (percentZoom < 0) ? -1 : 1;
    zoomPercentSetting += delta;
    calcScale("zoomByPercent");
  }

  private void setZoomParameters() {
    if (zoomPercentSetting < 5)
      zoomPercentSetting = 5;
    if (zoomPercentSetting > MAXIMUM_ZOOM_PERCENTAGE)
      zoomPercentSetting = MAXIMUM_ZOOM_PERCENTAGE;
    zoomPercent = (zoomEnabled) ? zoomPercentSetting : 100;
  }
  
  private void calcScale(String from) {
    setZoomParameters();
    scalePixelsPerAngstrom = scaleDefaultPixelsPerAngstrom * zoomPercent / 100;
  }
  

  void setZoomEnabled(boolean zoomEnabled) {
    if (this.zoomEnabled != zoomEnabled) {
      this.zoomEnabled = zoomEnabled;
      calcScale("setZoomEnabled");
    }
  }

  void setScaleAngstromsPerInch(float angstromsPerInch) {
    scalePixelsPerAngstrom = scaleDefaultPixelsPerAngstrom = 72 / angstromsPerInch;
  }

  /* ***************************************************************
   * PERSPECTIVE
   ****************************************************************/
  private boolean perspectiveDepth = true;
  private float cameraDepth = 3;
  private int cameraDistance = 1000; // prevent divide by zero on startup
  private float cameraDistanceFloat = 1000; // prevent divide by zero on startup

  void setPerspectiveDepth(boolean perspectiveDepth) {
    if (this.perspectiveDepth == perspectiveDepth)
      return;
    this.perspectiveDepth = perspectiveDepth;
    scaleFitToScreen();
  }

  boolean getPerspectiveDepth() {
    return perspectiveDepth;
  }

  /* ***************************************************************
   * SCREEN SCALING
   ****************************************************************/
  int width, height;
  int screenPixelCount;
  public float scalePixelsPerAngstrom;
  float scaleDefaultPixelsPerAngstrom;

  void setScreenDimension(int width, int height) {
    this.width = width;
    this.height = height;
  }

  private void setTranslationCenterToScreen() {
    // translate to the middle of the screen
    xFixedTranslation = width / 2;
    yFixedTranslation = height / 2;
    // 2005 02 22
    // switch to finding larger screen dimension
    // find smaller screen dimension
    screenPixelCount = width;
    boolean isByLarger = false;
    if (height != width && (isByLarger && height > width || !isByLarger && height < width))
      screenPixelCount = height;
    // ensure that rotations don't leave some atoms off the screen
    // note that this radius is to the furthest outside edge of an atom
    // given the current VDW radius setting. it is currently *not*
    // recalculated when the vdw radius settings are changed
    // leave a very small margin - only 1 on top and 1 on bottom
    if (screenPixelCount > 2)
      screenPixelCount -= 2;
    cameraScaleFactor = Float.MAX_VALUE;
  }
  
  private float defaultScaleToScreen(float radius) {
    /* 
     * 
     * the presumption here is that the rotation center is at pixel
     * (150,150) of a 300x300 window. rotationRadius is
     * a rough estimate of the furthest distance from the center of rotation
     * (but not including pmesh, special lines, planes, etc. -- just atoms)
     * 
     * also that we do not want it to be possible for the model to rotate
     * out of bounds of the applet. For internal spinning I had to turn
     * of any calculation that would change the rotation radius.  hansonr
     * 
     */
    return screenPixelCount / 2f / radius
        * getCameraScaleFactor(cameraDepth);
  }

  float cameraScaleFactor = Float.MAX_VALUE;
  private float getCameraScaleFactor(float depth) {
    if (!perspectiveDepth)
      return 1;
    /*
     *  Say you have a 300x300 applet. Then we really think of it as a
     *  300x300x300 cube and define the camera distance as some multiple
     *  (cameraDepth) of this 300 pixel "screen depth".
     *  If the camera is far, far away from the model, then there won't 
     *  be much perspective setting, and we don't need a scaling factor. 
     *  But if the camera is close in, then perspective will drive XY points
     *  near the camera outside the applet window unless we scale them down. 
     *  Note that the calculation below reduces to:
     *  
     *  scaleFactor = (cameraDepth + 0.5) / cameraDepth 
     *   
     *  or, simply
     *  
     *  scaleFactor = 1 + (0.5 / cameraDepth)
     *  
     *  I can find nothing anywhere in any code that sets cameraDepth other
     *  than its default value of 3, so I think scaleFactor is always 1.167
     *  prior to adding 0.02 "for luck".
     *    
     *  hansonr
     */

    if (cameraScaleFactor == Float.MAX_VALUE) {
      cameraDistance = (int) (depth * screenPixelCount);
      cameraDistanceFloat = cameraDistance;
      cameraScaleFactor = (cameraDistance + screenPixelCount / 2)
          / cameraDistanceFloat;
      // mth - for some reason, I can make the scaleFactor bigger in this
      // case. I do not know why, but there is extra space around the edges.
      // I have looked at it three times and still cannot figure it out
      // so just bump it up a bit.
      cameraScaleFactor += 0.02;
    }
    return cameraScaleFactor;
  }
    
  void scaleFitToScreen() {
    if (width == 0 || height == 0)
      return;
    setTranslationCenterToScreen();
    scaleDefaultPixelsPerAngstrom = defaultScaleToScreen(rotationRadius);
    calcScale("scaleFitToScreen rotrad=" + rotationRadius);
  }

  float perspectiveFactor(float z) {
    // all z's SHOULD be >= 0
    // so the more positive z is, the smaller the screen scale
    //new idea: phase out perspective depth when zoom is very large.
    //zoomPercent 1000 or larger starts removing this effect
    //we can go up to 200000
    //trouble with navitationMode is that it allows
    //z to go less than cameraDistance, thus producing
    //x-y EXPANSIONS rather than contractions. 
    
    float factor = (z <= 0? cameraDistanceFloat : cameraDistanceFloat / z);
    if (zoomPercent >= MAXIMUM_ZOOM_PERSPECTIVE_DEPTH)
      factor += (zoomPercent - MAXIMUM_ZOOM_PERSPECTIVE_DEPTH)/(MAXIMUM_ZOOM_PERCENTAGE - MAXIMUM_ZOOM_PERSPECTIVE_DEPTH) * (1 - factor);
    //System.out.println(z+" "+cameraDistanceFloat + " " + factor);
    return factor;  
  }
  
  short scaleToScreen(int z, int milliAngstroms) {
    if (milliAngstroms == 0)
      return 0;
    int pixelSize = (int) (milliAngstroms * scalePixelsPerAngstrom / 1000);
    if (perspectiveDepth)
      pixelSize *= perspectiveFactor(z);
    return (short) (pixelSize > 0 ? pixelSize : 1);
  }

  float scaleToPerspective(int z, float sizeAngstroms) {
    //DotsRenderer only
    return (perspectiveDepth ? sizeAngstroms * perspectiveFactor(z)
        : sizeAngstroms);
  }

  /* ***************************************************************
   * TRANSFORMATIONS
   ****************************************************************/

  final Matrix4f matrixTransform = new Matrix4f();
  private final Point3f point3fScreenTemp = new Point3f();
  private final Point3i point3iScreenTemp = new Point3i();
  private final Matrix4f matrixTemp = new Matrix4f();
  private final Vector3f vectorTemp = new Vector3f();

  /* ***************************************************************
   * RasMol has the +Y axis pointing down
   * And rotations about the y axis are left-handed
   * setting this flag makes Jmol mimic this behavior
   ****************************************************************/
  boolean axesOrientationRasmol = false;

  void setAxesOrientationRasmol(boolean axesOrientationRasmol) {
    this.axesOrientationRasmol = axesOrientationRasmol;
  }

  private final Vector3f perspectiveOffset = new Vector3f(0, 0, 0);

  synchronized void finalizeTransformParameters() {
    calcTransformMatrix();
    calcSlabAndDepthValues();
    // lock in the perspective so that when you change
    // centers there is no jump
    if(windowCentered) {
      matrixTransform.transform(rotationCenterDefault, pointT);
      matrixTransform.transform(fixedRotationCenter, pointT2);
      perspectiveOffset.sub(pointT, pointT2);
    }
    perspectiveOffset.x = xFixedTranslation; 
    perspectiveOffset.y = yFixedTranslation;
    if (windowCentered )
      perspectiveOffset.z = 0;

    /*
     * Note that the effect of this modification is restricted to the 
     * (undocumented) specialized circumstances when both
     * 
     * (a) set windowCentered is false (formerly "set frieda on")
     * 
     *  AND
     *
     * (b) the center has been changed to something other than the default
     * rotation center, either using "set picking center" followed by a user
     * click of an atom, or by a scripted "center (atom expression)".
     * 
     * This adjustment has no effect whatsoever on general use.
     * 
     * Bob Hanson 4/06
     *  
     */
  }
 
  synchronized private void calcTransformMatrix() {
    // you absolutely *must* watch the order of these operations
    matrixTransform.setIdentity();

    // first, translate the coordinates back to the center

    vectorTemp.set(fixedRotationCenter);
    matrixTemp.setZero();
    matrixTemp.setTranslation(vectorTemp);
    matrixTransform.sub(matrixTemp);
    // now, multiply by angular rotations
    // this is *not* the same as  matrixTransform.mul(matrixRotate);
    matrixTemp.set(matrixRotate);
    matrixTransform.mul(matrixTemp, matrixTransform);
    // we want all z coordinates >= 0, with larger coordinates further away
    // this is important for scaling, and is the way our zbuffer works
    // so first, translate an make all z coordinates negative
    vectorTemp.x = 0;
    vectorTemp.y = 0;
    vectorTemp.z = rotationRadius + cameraDistanceFloat
        / scalePixelsPerAngstrom;
    matrixTemp.setZero();
    matrixTemp.setTranslation(vectorTemp);
    if (axesOrientationRasmol)
      matrixTransform.add(matrixTemp); // make all z positive
    else
      matrixTransform.sub(matrixTemp); // make all z negative
    // now scale to screen coordinates
    matrixTemp.setZero();
    matrixTemp.set(scalePixelsPerAngstrom);
    if (!axesOrientationRasmol) {
      // negate y (for screen) and z (for zbuf)
      matrixTemp.m11 = matrixTemp.m22 = -scalePixelsPerAngstrom;
    }
    matrixTransform.mul(matrixTemp, matrixTransform);
    // note that the image is still centered at 0, 0 in the xy plane
    // all z coordinates are (should be) >= 0
    // translations come later (to deal with perspective)
  }

  Matrix4f getUnscaledTransformMatrix() {
    //for povray only
    Matrix4f unscaled = new Matrix4f();
    unscaled.setIdentity();
    vectorTemp.set(fixedRotationCenter);
    matrixTemp.setZero();
    matrixTemp.setTranslation(vectorTemp);
    unscaled.sub(matrixTemp);
    matrixTemp.set(matrixRotate);
    unscaled.mul(matrixTemp, unscaled);
    return unscaled;
  }

  void transformPoints(int count, Point3f[] angstroms, Point3i[] screens) {
    for (int i = count; --i >= 0;)
      screens[i].set(transformPoint(angstroms[i]));
  }

  void transformPoint(Point3f pointAngstroms, Point3i pointScreen) {
    pointScreen.set(transformPoint(pointAngstroms));
  }

  /** 
   * CAUTION! returns a POINTER TO A TEMPORARY VARIABLE
   * @param pointAngstroms
   * @return POINTER TO point3iScreenTemp
   */
  synchronized Point3i transformPoint(Point3f pointAngstroms) {
    matrixTransform.transform(pointAngstroms, point3fScreenTemp);
    return adjustedTemporaryScreenPoint();
  }
  
  Point3i adjustedTemporaryScreenPoint() {
    float z = (point3fScreenTemp.z - perspectiveOffset.z);
    if (z < cameraDistance) {
      if (Float.isNaN(point3fScreenTemp.z)) {
        //removed for extending pmesh to points and lines  BH 2/25/06 
        z = 1;
      } else if (z <= 0) {
        //just don't let z go past 1  BH 11/15/06
        z = 1;
      }
    }
    point3fScreenTemp.z = z;
    if (perspectiveDepth) {
      float perspectiveFactor = perspectiveFactor(z);
      point3fScreenTemp.x *= perspectiveFactor;
      point3fScreenTemp.y *= perspectiveFactor;
    }

    //higher resolution here for spin control. 

    point3fScreenTemp.x += perspectiveOffset.x;
    point3fScreenTemp.y += perspectiveOffset.y;

    point3iScreenTemp.x = (int) point3fScreenTemp.x;
    point3iScreenTemp.y = (int) point3fScreenTemp.y;
    point3iScreenTemp.z = (int) point3fScreenTemp.z;

    return point3iScreenTemp;
  }

  void unTransformPoint(Point3i screenPt, Point3f coordPt) {
    Point3f pt = new Point3f();
    pt.set(screenPt.x, screenPt.y, screenPt.z);
    pt.x -= perspectiveOffset.x;
    pt.y -= perspectiveOffset.y;
    if (perspectiveDepth) {
      float perspectiveFactor = perspectiveFactor(pt.z);
      pt.x /= perspectiveFactor;
      pt.y /= perspectiveFactor;
    }
    pt.z += perspectiveOffset.z;
    Matrix4f m = new Matrix4f();
    m.invert(matrixTransform);
    m.transform(pt, coordPt);
  }

  void transformVector(Vector3f vectorAngstroms, Vector3f vectorTransformed) {
    matrixTransform.transform(vectorAngstroms, vectorTransformed);
  }

  /////////// rotation center ////////////
  
  //from Frame:
  
  boolean windowCentered;
  
  boolean isWindowCentered() {
    return windowCentered;
  }

  void setWindowCentered(boolean TF) {
    windowCentered = TF;
  }

  /* This fixes the growing model clipping bug.
   * Best guess: the rotation radius does not update as the model grows, and
   * scaling eventually fails to keep the model within the slab.
   */
  void setDefaultRotation() {
    rotationCenterDefault = g3dsys.getBoundingBoxCenter();
    setFixedRotationCenter(rotationCenterDefault);
    rotationRadius = rotationRadiusDefault = g3dsys.calcRotationRadius(rotationCenterDefault);
    windowCentered = true;
  }

  Point3f getRotationCenter() {
    return fixedRotationCenter;
  }

  float getRotationRadius() {
    return rotationRadius;
  }

  private void setRotationCenterAndRadiusXYZ(Point3f newCenterOfRotation,
                                             boolean andRadius) {
    if (newCenterOfRotation == null) {
      setFixedRotationCenter(rotationCenterDefault);
      rotationRadius = rotationRadiusDefault;
      return;
    }
    setFixedRotationCenter(newCenterOfRotation);
    if (andRadius && windowCentered)
      rotationRadius = g3dsys.calcRotationRadius(fixedRotationCenter);
  }

  /* ***************************************************************
   * SLAB
   ****************************************************************/

  /*
   slab is a term defined and used in rasmol.
   it is a z-axis clipping plane. only atoms behind the slab get rendered.
   100% means:
   - the slab is set to z==0
   - 100% of the molecule will be shown
   50% means:
   - the slab is set to the center of rotation of the molecule
   - only the atoms behind the center of rotation are shown
   0% means:
   - the slab is set behind the molecule
   - 0% (nothing, nada, nil, null) gets shown
   */


  boolean slabEnabled = true;

  int slabPercentSetting = 100;
  int depthPercentSetting = 0;

  public int slabValue;
  public int depthValue;

  int getSlabPercentSetting() {
    return slabPercentSetting;
  }

  /**
   * Gets depth percent
   * @return returns current depth percent
   */
  // added for symmetry 2/15/2007 --C.T.
  int getDepthPercentSetting() {
    return depthPercentSetting;
  }
  
  /**
   * Increments slab percentage
   * @param percentage the amount to increase by
   */
  void slabByPercentagePoints(int percentage) {
    slabPercentSetting += percentage;
    if (slabPercentSetting < 1)
      slabPercentSetting = 1;
    else if (slabPercentSetting > 100)
      slabPercentSetting = 100;
    if (depthPercentSetting >= slabPercentSetting)
      depthPercentSetting = slabPercentSetting - 1;
  }

  /**
   * Increments depth percentage
   * @param percentage the amount to increase by
   */
  void depthByPercentagePoints(int percentage) {
    depthPercentSetting += percentage;
    if (depthPercentSetting < 0)
      depthPercentSetting = 0;
    else if (depthPercentSetting > 99)
      depthPercentSetting = 99;
    if (slabPercentSetting <= depthPercentSetting)
      slabPercentSetting = depthPercentSetting + 1;
  }

  /**
   * Increments slab/depth percentage at once
   * @param percentage the amount to increase by
   */
  void slabDepthByPercentagePoints(int percentage) {
    if (percentage > 0) {
      if (slabPercentSetting + percentage > 100)
        percentage = 100 - slabPercentSetting;
    } else {
      if (depthPercentSetting + percentage < 0)
        percentage = 0 - depthPercentSetting;
    }
    slabPercentSetting += percentage;
    depthPercentSetting += percentage;
  }

  /**
   * Sets a specific slab percentage
   * @param percentSlab
   */
  void slabToPercent(int percentSlab) {
    slabPercentSetting = percentSlab < 1 ? 1 : percentSlab > 100 ? 100
        : percentSlab;
    if (depthPercentSetting >= slabPercentSetting)
      depthPercentSetting = slabPercentSetting - 1;
  }

  /**
   * Enable/disable slab and depth
   * @param slabEnabled
   */
  public void setSlabEnabled(boolean slabEnabled) {
    this.slabEnabled = slabEnabled;
  }


  /**
   * Sets a specific depth percentage
   * @param percentDepth
   */
  // depth is an extension added by OpenRasMol
  // it represents the 'back' of the slab plane
  void depthToPercent(int percentDepth) {
    depthPercentSetting = percentDepth < 0 ? 0 : percentDepth > 99 ? 99
        : percentDepth;
    if (slabPercentSetting <= depthPercentSetting)
      slabPercentSetting = depthPercentSetting + 1;
  }

  /**
   * Converts slab and depth percentages into the actual (zbuffer?) values
   * used by g3d. 
   */
  //XXX Note that the g3d documentation lies about using percentages
  private void calcSlabAndDepthValues() {
    slabValue = 0;
    depthValue = Integer.MAX_VALUE;
    if (slabEnabled) {
      // miguel 24 sep 2004 -- the comment below does not seem right to me
      // I don't think that all transformed z coordinates are negative
      // any more
      //
      // all transformed z coordinates are negative
      // a slab percentage of 100 should map to zero
      // a slab percentage of 0 should map to -diameter
      int radius = (int) (rotationRadius * scalePixelsPerAngstrom);
      slabValue = ((100 - slabPercentSetting) * 2 * radius / 100)
          + cameraDistance;
      depthValue = ((100 - depthPercentSetting) * 2 * radius / 100)
          + cameraDistance;
    }
  }
  
}
