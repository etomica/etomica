package g3dsys.control;

import g3dsys.images.Figure;

import javax.vecmath.Point3f;

/**
 *	Class that stores figures and delegates draw commands to them 
 */

class FigureManager {
	
	private java.util.Collection figs; //holds Figures to be drawn.
	
	private long idCount; //for giving Figures IDs when they are added
	
	/* Store molSpace size in Angstroms; needed for proper scaling
	 * These are determined by the contents of the model, and change as atoms
	 * are added and deleted.
	 */
	private Point3f min = new Point3f(0,0,0);
	private Point3f max = new Point3f(0,0,0);
	
	private G3DSys gsys;
	
	public FigureManager(G3DSys g) {
		figs = new java.util.ArrayList();
		idCount = 0;
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
		gsys.setCenterOfRotation(new javax.vecmath.Point3f(
				(min.x+max.x)/2, (min.y+max.y)/2, (min.z+max.z)/2));
		
		long start = System.currentTimeMillis();
		for(java.util.Iterator iter = figs.iterator(); iter.hasNext();) {
			Figure f = (Figure)iter.next();
			if( f != null ) f.draw();
		}
		System.out.println("render time "+(System.currentTimeMillis()-start));
		/*
		 * The same on this level: rendering time increases as window size grows,
		 * but after clipping it repaints faster despite no change in render
		 * time
		 */
	}
	
	public long addFig(Figure f) {
		long id = addFigNoRescale(f);
		gsys.recalcPPA();
		return id;
	}

	/**
	 * Stores an additional figure, expanding model bounds as needed
	 * @param f the Figure to add
	 */
	public long addFigNoRescale(Figure f) {
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

		f.setID(idCount);
		figs.add(f);
		
		return idCount++;
	}

	public Figure removeFig(long id) {
		Figure f = removeFigNoRescale(id);
		shrinkModel();
		gsys.recalcPPA();
		return f;
	}
	
	/**
	 * Remove a Figure from the system without resizing.
	 * Useful when doing batch removals from a large system, but the user
	 * must manually call shrinkModel at the end.
	 * @param id the ID of the figure to remove
	 * @return the removed figure (or null)
	 */
	public Figure removeFigNoRescale(long id) {
		for(java.util.Iterator iter = figs.iterator(); iter.hasNext();) {
			Figure f = (Figure)iter.next();
			if(f.getID() == id) {
				figs.remove(f);
				return f;
			}
		}
		return null;
	}

	/**
	 * Resizes min and max Angstrom values for the model based on what Figures
	 * are still present
	 * Uses iterators. Removing n figures in succession can be quadratic.
	 */
	public void shrinkModel() {
		float xn,yn,zn; //miN
		float xm,ym,zm; //Max
		java.util.Iterator i = figs.iterator();
		if(!i.hasNext()) return;

		Figure f = (Figure)i.next();
		//System.out.println(""+f.getID()+": "+f.getX()+","+f.getY()+","+f.getZ());
		float offset = f.getD()/2;
		xn = f.getX() - offset; yn = f.getY() - offset; zn = f.getZ() - offset;
		xm = f.getX() + offset; ym = f.getY() + offset; zm = f.getZ() + offset;
		
		while(i.hasNext()) {
			f = (Figure)i.next();
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

	/**
	 * Gets an array of longs for the current collection of Figures.
	 * User must refresh stale data if model is modified afterwards.
	 * @return the array of current Figure IDs.
	 */
	public long[] getFigs() {
		Figure[] figarr = (Figure[])figs.toArray();
		long[] idarr = new long[figarr.length];
		for(int i=0; i<figarr.length; i++)
			idarr[i] = figarr[i].getID();
		return idarr;
	}	
	
}
