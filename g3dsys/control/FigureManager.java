package g3dsys.control;

import g3dsys.images.*;

import javax.vecmath.Point3f;
import javax.vecmath.Tuple3f;

/**
 *	Class that stores figures and delegates draw commands to them 
 */

class FigureManager {
	
    private Figure[] figs;
	
    private int idMax; //for giving Figures IDs when they are added
	
    /* Store molSpace size in Angstroms; needed for proper scaling
     * These are determined by the contents of the model, and change as atoms
     * are added and deleted.
     */
    private Point3f min = new Point3f(0,0,0);
    private Point3f max = new Point3f(0,0,0);
    private Point3f tempP = new Point3f();

    private G3DSys gsys;

    public FigureManager(G3DSys g) {
        //figs = new java.util.HashSet();
        figs = new Figure[0];
        idMax = -1;
        gsys = g;
        images = new ImageShell(gsys);
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
	public void draw() {
		// persistent center of rotation in the middle of the model
//        tempP.x = (min.x+max.x)/2;
//        tempP.y = (min.y+max.y)/2;
//        tempP.z = (min.z+max.z)/2;
//        gsys.setCenterOfRotation(tempP);
		
        for (int j=0; j<idMax+1; j++) {
            figs[j].draw();
		}
	}
	
	/**
	 * Stores an additional figure, expanding model bounds as needed
	 * @param f the Figure to add
	 */
    public void addFig(Figure f) {
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

	/**
	 * Remove a Figure from the system without resizing.
	 * Useful when doing batch removals from a large system, but the user
	 * must manually call shrinkModel at the end.
	 * @param id the ID of the figure to remove
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
        gsys.recalcPPA();
    }

    /**
     * Finds the furthest distance in the model; in our dynamic model from
     * one corner to the other. For rotation, assumed about the origin.
     * This will always be a safe distance, regardless of the rotation point.
     * @param center the reference point; currently ignored
     * @return the furthest distance found, in Angstroms
     */
    public float calcRotationRadius(Point3f center) {
      //real radius is /2, but still has clipping when atoms are near the boundary
      //and /1 shrinks the model too much at 100% zoom; /1.5f compromise
      //note that if an atom radius causes it to extend well outside the box
      //there may still be clipping
      if(imagesOn) {
        //return a different radius to account for image shell
        //getLayers*2 for symmetry, +1 for original center image
        return (min.distance(max)/1.5f)*(2*images.getLayers()+1);
      }
      else return min.distance(max)/1.5f;
    }

    public Point3f getBoundingBoxCenter() {
      return new Point3f(0,0,0);
    }

    public Point3f getAverageAtomPoint() {
      Point3f average = new Point3f(0,0,0);
      for (int i = figs.length; --i >= 0;)
        average.add(figs[i].getPoint()); //nulls in figs?
      average.scale(1f / figs.length);
      return average;
    }

    /* **********************************************************
     * Image shell code 
     ************************************************************/

    ImageShell images;
    private boolean imagesOn = false;
    
    /**
     * For use only by ImageShell class for iteration
     * @return returns the array of current Figures
     */
    public Figure[] getFigs() {
      return figs;
    }

    //stores old values while image shell is on; may cause problems
    //with dynamic system...
    Point3f oldmin = new Point3f();
    Point3f oldmax = new Point3f();
    
    public boolean isEnableImages() { return imagesOn; }
    public void setEnableImages(boolean b) {
      //save some array overhead with these checks
      if(imagesOn && b) return; //no change
      if(!imagesOn && !b) return; //no change
      if(!imagesOn && b) addFig(images); //toggle on
      if(imagesOn && !b) removeFig(images); //toggle off
      imagesOn = b;
    }

}
