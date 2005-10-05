package etomica.graphics2;


import etomica.atom.Atom;
import etomica.atom.AtomFilter;
import etomica.atom.AtomFilterStatic;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomTypeWell;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polytope;
import etomica.phase.Phase;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * Converts all drawable information in a Phase object into graphics primitives
 *  using the Renderable interface. 
 * It does not handle mouse or keyboard - it's up to the implamentation to do it. 
 */
public final class SceneManager implements  java.io.Serializable
{
    
    public SceneManager() 
    {
        setAtomFilter(AtomFilterStatic.ACCEPT_ALL);
        setScale(1.0);
 
    	colorScheme = new ColorSchemeByType();
    	
   }

    public void setRenderer( Renderable r )
    {
    	renderer = r;
    }

   
    public void updateAtomPositions()
    {
        if(!isVisible() || getPhase() == null) return;
	    if(!tablesInitialized) 
		  	  initInternalTables();
    	for ( int i=0; i<sphereCores.length; i++ )
    	{
   		      Atom a = sphereCores[i];
   		      int c = colorScheme.atomColor(a);
   		      Vector r = a.coord.position();

   		      Renderable.Shape shp = core_shapes[i];
   		      shp.setPosition( r );
   		      shp.setColor( c );
		}
    }
    
    public void setPhase(Phase phase) {
    	this.phase = phase;
    }
    
    public Phase getPhase() {
    	return phase;
    }
    
    public boolean isVisible() {
    	return true;
    }
    
	public ColorScheme getColorScheme() {
		return colorScheme;
	}
	public void setColorScheme(ColorScheme colorScheme) {
		this.colorScheme = colorScheme;
	}
	
    public void setAtomFilter(AtomFilter filter) {atomFilter = filter;} 

    public void setWriteScale(boolean s) {writeScale = s;}
    public boolean getWriteScale() {return(writeScale);}

    public void setQuality(int q) {
      if(q > DRAW_QUALITY_VERY_HIGH)
        q -= DRAW_QUALITY_MAX;
      if(q < DRAW_QUALITY_VERY_LOW)
        q += DRAW_QUALITY_MAX;
      quality = q;
    }
    public int getQuality() {return(quality);}
    
    public void setDrawBoundary(int b) {
      if(b>DRAW_BOUNDARY_ALL)
        b-=DRAW_BOUNDARY_MAX;
      else if(b<DRAW_BOUNDARY_NONE)
        b+=DRAW_BOUNDARY_MAX;
      drawBoundary = b;
    }
    public int getDrawBoundary() {return drawBoundary;}

    
    public synchronized void initInternalTables() 
    {
    	try
		{
    		int countSphereCores = 0, countSphereWells = 0, countSphereRotators = 0;
    		int countWalls = 0, countAll = 0;
    		float vertAll[];
    		Atom atoms[];
    		
    		countAll = phase.getSpeciesMaster().node.leafAtomCount();
    		
    		if(countAll==0) return;
    		
    		vertAll = new float[countAll*3];
    		atoms = new Atom[countAll];
    		
    		int i = 0;
    		AtomIteratorLeafAtoms iter = new AtomIteratorLeafAtoms(phase);
    		iter.reset();
    		while(iter.hasNext()) {
    			Atom a = iter.nextAtom();
    			atoms[i/3] = a;
    			vertAll[i] = (float)a.coord.position().x(0);// + drawExpansionShiftX;
    			vertAll[i+1] = (float)a.coord.position().x(1);// + drawExpansionShiftY;
    			vertAll[i+2] = (float)a.coord.position().x(2);// + drawExpansionShiftZ;
    			if(a.type instanceof AtomTypeOrientedSphere) countSphereRotators++;
    			if(a.type instanceof AtomTypeWell) countSphereWells++;
    			if(a.type instanceof AtomTypeSphere) countSphereCores++;
    			i += 3;
    		}
    		
    		sphereCores = new Atom[countSphereCores];
    		sphereWells = new Atom[countSphereWells];
    		sphereRotators = new Atom[countSphereRotators];
    		walls = new Atom[countWalls];
    		vertSphereCores = new float[countSphereCores*3];
    		vertSphereWellBase = new int[countSphereWells];
    		vertSphereRotatorBase = new int[countSphereRotators];
    		vertWalls = new float[countWalls*3];
    		for(int j=0,k=0,l=0,m=0; (j/3) < atoms.length; j+=3) {
    			if(atoms[j/3].type instanceof AtomTypeSphere) {
    				sphereCores[m/3] = atoms[j/3];
    				vertSphereCores[m] = vertAll[j];
    				vertSphereCores[m+1] = vertAll[j+1];
    				vertSphereCores[m+2] = vertAll[j+2];
    				m+=3;
    			}
    			if(atoms[j/3].type instanceof AtomTypeOrientedSphere) {
    				sphereRotators[k] = atoms[j/3];
    				vertSphereRotatorBase[k] = m-3;
    				k++;
    			}
    			if(atoms[j/3].type instanceof AtomTypeWell) {
    				sphereWells[l] = atoms[j/3];
    				vertSphereWellBase[l] = m-3;
    				l++;
    			}
    		}
    		
    		// Create graphical object for sphere cores
    		core_shapes = new Renderable.Shape[ countSphereCores ];
    		Vector3D scale = new Vector3D();
    		for ( int j=0; j<countSphereCores; j++ )
    		{
    			// Create object (it's a unit sphere)
    			Renderable.Sphere shp = renderer.createSphere();
    			core_shapes[j] = shp;
    			
    			// Scale to atom's given size
    			Atom a = sphereCores[j];
    			double radius = ((AtomTypeSphere)a.type).diameter( a );
    			scale.E( radius );
    			shp.setScale( scale );
    		}
    		
    		// Create graphical object for the boundary
    		Boundary bnd = phase.boundary();
    		Polytope shape = bnd.getShape();
    	    LineSegment[] edges = shape.getEdges();
    	    Vector shift = phase.space().makeVector();
    	    shift.Ea1Tv1( +0.5, bnd.dimensions() );
    	    
    	    Renderable.Polyline poly = renderer.createPoly();
	    	Vector from = phase.space().makeVector();
	    	Vector to = phase.space().makeVector();
    	    for(i=0; i<edges.length; i++) 
    	    {  	        
    	    	from.Ev1Pv2( edges[i].getVertices()[0], shift );
    	    	to.Ev1Pv2( edges[i].getVertices()[1], shift );
    	    	poly.appendLine( from, to );    	        
    	    }
		}
    	finally {
    		tablesInitialized = true;
    	}
    }
    
	public double getScale() {
		return scale;
	}
	public void setScale(double scale) {
		this.scale = scale;
	}
	
	public Atom[] getSelectedAtoms() {
		return selectedAtoms;
	}
	public void setSelectedAtoms(Atom[] selectedAtoms) {
		this.selectedAtoms = selectedAtoms;
	}
	
       
    static final int DRAW_QUALITY_VERY_LOW = 0;
    static final int DRAW_QUALITY_LOW = 1;
    static final int DRAW_QUALITY_NORMAL = 2;
    static final int DRAW_QUALITY_HIGH = 3;
    static final int DRAW_QUALITY_VERY_HIGH = 4;
    static final int DRAW_QUALITY_MAX = 5;
    //Boundary Constants
    static final int DRAW_BOUNDARY_NONE = 0;
    static final int DRAW_BOUNDARY_OUTLINE = 1;
    static final int DRAW_BOUNDARY_SHELL = 2;
    static final int DRAW_BOUNDARY_ALL = 3;
    static final int DRAW_BOUNDARY_MAX = 4;
    

    protected Renderable.Shape[] core_shapes;
    
    protected final AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms();
    protected AtomFilter atomFilter;
    protected double scale;
    protected Phase phase;
    
    protected ColorScheme colorScheme;
    protected Atom[] selectedAtoms = new Atom[1];
    protected Renderable renderer;
    private boolean tablesInitialized = false;

//  The groups of atoms
    private Atom sphereCores[];
    private Atom sphereWells[];
    private Atom sphereRotators[];
    private Atom walls[];
//    The verticies of said atoms
    private float vertSphereCores[];
    private int vertSphereWellBase[];
    private int vertSphereRotatorBase[];
    private float vertWalls[];


   /**
    * Variable specifying whether a line tracing the boundary of the display should be drawn
    * Default value is <code>BOUNDARY_OUTLINE</code>
    */
    int drawBoundary = DRAW_BOUNDARY_OUTLINE;

    /**
     * Variable that sets the quality of the rendered image.
     */
    int quality = DRAW_QUALITY_NORMAL;

    /** 
     * Flag to indicate if value of scale should be superimposed on image
     */
    boolean writeScale = false;
    
    /**
     *  Sets the quality of the rendered image, false = low, true = high
      */
    boolean highQuality = false;

} //end of ConfigurationCanvas class

