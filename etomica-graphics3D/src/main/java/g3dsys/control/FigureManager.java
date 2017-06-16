/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package g3dsys.control;

import g3dsys.images.Ball;
import g3dsys.images.Bond;
import g3dsys.images.Figure;
import g3dsys.images.ImageShell;
import org.jmol.util.Point3f;

/**
 *	Class that stores figures and delegates draw commands to them 
 */

class FigureManager {

  ImageShell images;
  //stores old values while image shell is on; may cause problems
  //with dynamic system...
  Point3f oldmin = new Point3f();
  Point3f oldmax = new Point3f();
  private Figure[] figs;
  private int idMax; //for giving Figures IDs when they are added
  /* Store molSpace size in Angstroms; needed for proper scaling
   * These are determined by the contents of the model, and change as atoms
   * are added and deleted.
   */
  private Point3f min = Point3f.new3(0,0,0);
  private Point3f max = Point3f.new3(0,0,0);
  private G3DSys gsys;
  private boolean imagesOn = false;
  private boolean wireframe = false;
  public FigureManager(G3DSys g) {
    //figs = new java.util.HashSet();
    figs = new Figure[0];
    idMax = -1;
    gsys = g;
    images = new ImageShell(gsys);
    addFig(images); //imageshell is the first figure in the array
  }

  /** Get the depth of the model in Angstroms
   *  @return depth of the model in Angstroms */
  public float getDepth() { return max.z - min.z; }

  /** Get the width of the model in Angstroms
   * @return width of the model in Angstroms */
  public float getWidth() { return max.x - min.x; }

  /** Get the height of the model in Angstroms
   * @return height of the model in Angstroms */
  public float getHeight() { return max.y - min.y; }

  /** Get the model's minimum x value in Angstroms
   * @return minimum x value in Angstroms */
  public float getMinX() { return min.x; }

  /** Get the model's minimum y value in Angstroms
   * @return minimum y value in Angstroms */
  public float getMinY() { return min.y; }

  /** Get the model's minimum z value in Angstroms
   * @return minimum z value in Angstroms */
  public float getMinZ() { return min.z; }

  /** Get the model's maximum x value in Angstroms
   * @return maximum x value in Angstroms */
  public float getMaxX() { return max.x; }

  /** Get the model's maximum y value in Angstroms
   * @return maximum y value in Angstroms */
  public float getMaxY() { return max.y; }

  /** Get the model's maximum z value in Angstroms
   * @return maximum z value in Angstroms */
  public float getMaxZ() { return max.z; }

  /** Dispatches draw commands to all stored Figures */
  public synchronized void draw() {
    for (int j=0; j<idMax+1; j++) {
      if(figs[j] instanceof Ball && wireframe) continue;
      figs[j].draw();
    }
  }

  /**
   * Stores an additional figure, expanding model bounds as needed
   * @param f the Figure to add
   */
  public synchronized void addFig(Figure f) {
    if (f.getID() > -1) {
      throw new IllegalArgumentException("figure is already here");
    }

    f.setID(++idMax);
    if (figs.length < idMax+1) {
      // no room in the array.  reallocate the array with an extra cushion.
      Figure[] newFigsArray = new Figure[(int)(idMax*1.5)+50];
      System.arraycopy(figs, 0, newFigsArray, 0, figs.length);
      figs = newFigsArray;
    }
    figs[idMax] = f;
  }

  /* **********************************************************
   * Image shell code 
   ************************************************************/

  /**
   * Remove a Figure from the system without resizing.
   * Useful when doing batch removals from a large system, but the user
   * must manually call shrinkModel at the end.
   * @param f the figure to remove
   * @return the removed figure (or null)
   */
  public synchronized Figure removeFig(Figure f) {
    if (f.getID() > idMax || figs[f.getID()] != f) {
      throw new IllegalArgumentException("Don't know about "+f);
    }
    int oldID = f.getID();
    if (oldID < idMax && idMax > 0) {
      figs[oldID] = figs[idMax];
      figs[oldID].setID(oldID);
      figs[idMax] = null;
    }
    else {
      figs[oldID] = null;
    }
    idMax--;
    if (idMax > 100 && idMax < figs.length/2) {
      Figure[] newFigsArray = new Figure[idMax+50];
      System.arraycopy(figs, 0, newFigsArray, 0, idMax+1);
      figs = newFigsArray;
    }
    f.setID(-1);
    return f;
  }

  public void setBoundingBox(float minx, float miny, float minz,
      float maxx, float maxy, float maxz) {
    min.x = minx; min.y = miny; min.z = minz;
    max.x = maxx; max.y = maxy; max.z = maxz;
  }

  /**
   * Finds the furthest distance in the model; in our dynamic model from
   * one corner to the other. For rotation, assumed about the origin.
   * This will always be a safe distance, regardless of the rotation point.
   * @param center the reference point; currently ignored
   * @return the furthest distance found, in Angstroms
   */
  public float calcRotationRadius(Point3f center) {
    /* Real radius is /2, but still has clipping when atoms are near the
     * boundary, and /1 shrinks the model too much at 100% zoom; /1.5f
     * compromise. Note that if an atom radius causes it to extend well
     * outside the box there may still be clipping.
     */
    float d = min.distance(max);
    if (d == 0) {
      d = 10;
    }
    if(imagesOn) {
      //return a different radius to account for image shell
      //getLayers*2 for symmetry, +1 for original center image
      return (d/1.5f)*(2*images.getLayers()+1);
    }
    return d/1.5f;
  }

  public Point3f getBoundingBoxCenter() {
    return new Point3f();
  }

  public synchronized Point3f getAverageAtomPoint() {
    Point3f average = Point3f.new3(0,0,0);
    int ballCount = 0;
    for (int i = figs.length; --i >= 0;)
      if (figs[i] instanceof Ball) {
        ballCount++;
        average.add(((Ball)figs[i]).getPoint()); //nulls in figs?
      }
    average.scale(1f / ballCount);
    return average;
  }

  /**
   * For use only by ImageShell class for iteration
   * @return returns the array of current Figures
   */
  public Figure[] getFigs() {
    return figs;
  }

  public boolean isEnableImages() { return imagesOn; }

  public void setEnableImages(boolean b) {
    //save some array overhead with these checks
    if(imagesOn && b) return; //no change
    if(!imagesOn && !b) return; //no change
    if(!imagesOn && b) images.setDrawable(true); //toggle on
    if(imagesOn && !b) images.setDrawable(false); //toggle off
    imagesOn = b;
    gsys.setDefaultRotation();
    gsys.fastRefresh();
  }

  /**
   * Gets the boundary vectors for the phase; used by ImageShell
   * @return returns an array representation of the vectors
   */
  public double[] getBoundaryVectors() {
    return images.getBoundaryVectors();
  }

  /**
   * Sets the boundary vectors for the phase; used by ImageShell
   *
   * @param values an array representation of the boundary vectors
   */
  public void setBoundaryVectors(double[] values) {
    images.setBoundaryVectors(values);
  }

  /**
   * Gets the number of image shell layers
   *
   * @return returns the number of layers
   */
  public int getLayers() {
    return images.getLayers();
  }

  /**
   * Sets the number of image shell layers
   * @param n the number of layers
   */
  public void setLayers(int n) {
    images.setLayers(n);
    gsys.setDefaultRotation();
    gsys.fastRefresh();
  }

  /**
   * Gets the boundary drawing style
   * @return returns the boundary drawing style
   */
  public int getDrawBoundaryType() { return images.getDrawBoundaryType();
  }

  /**
   * Sets the boundary drawing style
   *
   * @param b the boundary drawing style to use
   */
  public void setDrawBoundaryType(int b) {
    images.setDrawBoundaryType(b); }

  /**
   * Cycles to the next boundary drawing style
   */
  public void cycleDrawBoundaryType() { images.cycleDrawBoundaryType(); }

  /**
   * Toggle wireframe on/off.
   * Turn atom drawing off, bond types to wireframe; or
   * turn atom drawing on, bond types to cylinder.
   */
  public synchronized void toggleWireframe() {
    wireframe = !wireframe;
    for(int i=0; i<figs.length; i++) {
      if(figs[i] instanceof Bond) ((Bond)figs[i]).setBondType(
          (wireframe ? Bond.WIREFRAME : Bond.CYLINDER));
    }
    images.setWireFrame(wireframe);
  }

  public IndexIterator getBoundaryVectorsIterator() {
    return images.getBoundaryVectorsIterator();
  }

  public void setBoundaryVectorsIterator(IndexIterator i) {
    images.setBoundaryVectorsIterator(i);
  }

}
