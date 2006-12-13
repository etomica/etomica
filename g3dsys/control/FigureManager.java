package g3dsys.control;

import g3dsys.images.Figure;

import javax.vecmath.Point3f;

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
        tempP.x = (min.x+max.x)/2;
        tempP.y = (min.y+max.y)/2;
        tempP.z = (min.z+max.z)/2;
        gsys.setCenterOfRotation(tempP);
		
        for (int j=0; j<idMax+1; j++) {
            figs[j].draw();
		}
	}
	
	public void addFig(Figure f) {
		addFigNoRescale(f);
		shrinkModel();
		gsys.recalcPPA();
	}

	/**
	 * Stores an additional figure, expanding model bounds as needed
	 * @param f the Figure to add
	 */
	public void addFigNoRescale(Figure f) {
        if (f.getID() > -1) {
            throw new IllegalArgumentException("figure is already here");
        }

        //possible new min and max values due to the new figure
		float curMinX = f.getX() - f.getD()/2;
		float curMaxX = f.getX() + f.getD()/2;
		float curMinY = f.getY() - f.getD()/2;
		float curMaxY = f.getY() + f.getD()/2;
		float curMinZ = f.getZ() - f.getD()/2;
		float curMaxZ = f.getZ() + f.getD()/2;
		
		// if new coordinates exceed the old, expand space and store old max
		if(curMinX < min.x) { min.x = curMinX; }
		if(curMaxX > max.x) { max.x = curMaxX; }
		if(curMinY < min.y) { min.y = curMinY; }
		if(curMaxY > max.y) { max.y = curMaxY; }
		if(curMinZ < min.z) { min.z = curMinZ; }
		if(curMaxZ > max.z) { max.z = curMaxZ; }

        f.setID(++idMax);
        if (figs.length < idMax+1) {
          // no room in the array.  reallocate the array with an extra cushion.
          Figure[] newFigsArray = new Figure[(int)(idMax*1.5)+50];
          System.arraycopy(figs, 0, newFigsArray, 0, figs.length);
          figs = newFigsArray;
        }
        figs[idMax] = f;
	}

	public Figure removeFig(Figure f) {
		Figure r = removeFigNoRescale(f);
		shrinkModel();
		gsys.recalcPPA();
		return r;
	}
	
	/**
	 * Remove a Figure from the system without resizing.
	 * Useful when doing batch removals from a large system, but the user
	 * must manually call shrinkModel at the end.
	 * @param id the ID of the figure to remove
	 * @return the removed figure (or null)
	 */
	public Figure removeFigNoRescale(Figure f) {
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

	/**
	 * Resizes min and max Angstrom values for the model based on what Figures
	 * are still present
	 * Uses iterators. Removing n figures in succession can be quadratic.
	 */
	public void shrinkModel() {
		float xn,yn,zn; //miN
		float xm,ym,zm; //Max
		if(idMax < 0) return;

		Figure f = figs[0];
		//System.out.println(""+f.getID()+": "+f.getX()+","+f.getY()+","+f.getZ());
		float offset = f.getD()/2;
		xn = f.getX() - offset; yn = f.getY() - offset; zn = f.getZ() - offset;
		xm = f.getX() + offset; ym = f.getY() + offset; zm = f.getZ() + offset;
		
		for (int i=1; i<idMax+1; i++) {
			f = figs[i];
			offset = f.getD()/2;
			//System.out.println(""+f.getID()+": "+f.getX()+","+f.getY()+","+f.getZ());
			if( f.getX()-offset < xn ) xn = f.getX()-offset;
			if( f.getX()+offset > xm ) xm = f.getX()+offset;
			if( f.getY()-offset < yn ) yn = f.getY()-offset;
			if( f.getY()+offset > ym ) ym = f.getY()+offset;
			if( f.getZ()-offset < zn ) zn = f.getZ()-offset;
			if( f.getZ()+offset > zm ) zm = f.getZ()+offset;
		}
		
		min.x = xn; min.y = yn; min.z = zn;
		max.x = xm; max.y = ym; max.z = zm;
	}
}
