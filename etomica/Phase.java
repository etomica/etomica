package simulate;
import java.awt.*;//Because this class subclasses Component.
import java.io.*;//for Serializable and Stream.
import java.util.*; //for type Vector
import java.beans.*;//for PropertyChangeSupport and other beans convenience classes.
import java.util.*; //for neighborList
import java.awt.event.*;

/**
 * A Phase comprises a Space, one or more Species, and Potentials that characterize
 * how the species interact.  It may be viewed as a data structure that organizes
 * this information and which may be acted upon by an Integrator (which is the object
 * that contains information about how to do dynamics and ensemble sampling).
 * The Phase is responsible for knowing how to report instantaneous values (as opposed
 * to averages) of various mechanical and structural quantities (such as its energy
 * or a histogram of atom distances), and it knows how to draw itself to the screen.
 *
 * A Phase is a bean.  It is placed inside a Simulation, and it contains one Space, 
 * one or more Species, and one or more Potential1 and Potential2 objects.
 *
 * @author David Kofke
 * @author C. Daniel Barnes
 *
 * @see Simulation
 * @see Space
 * @see Species
 * @see Potential1
 * @see Potential2
 */
public final class Phase extends Container {

 /**
  * Conversion factor from simulation unit of length (Angstroms) to screen pixels (default = 300).
  * It is a class variable, so any instance of Phase that changes its value will
  * affect how all Phases are drawn to the screen.  The instance variable maySetTO_PIXELS
  * determines whether a phase has permission to change TO_PIXELS.
  */
  public static double TO_PIXELS = 300.;
 
 /**
  * Flag specifying whether phase computes energies, etc., using neighbor lists.  
  * Neighbor list feature is under revision, so presently their use should be avoided.
  */
  public  boolean useNeighborList;
 
 /**
  * Flag specifying whether a line tracing the boundary of the phase should be 
  * included when drawing the phase to the screen.
  * Default value is <code>false</code>
  */
  private boolean drawBoundingBox = false;
  
 /**
  * Flag indicating whether this phase has permission to change the value of TO_PIXELS
  * Default value is <code>true</code>.
  *
  * @see #TO_PIXELS
  */
  private boolean maySetTO_PIXELS = true;
  
 /**
  * Number of periodic-image shells to be drawn when drawing this phase to the
  * screen.  Default value is 0.
  *
  * @see #paint
  */
  private int imageShells = 0;
 
 /**
  * The nominal scaling factor that determines the size of this phase when drawn to the screen.
  * 
  * @see #paint
  */
  private double nominalScale = 1.0;
  
 /**
  * Total number of species contained in this phase.
  */
  int speciesCount=0;
  private transient int[] origin;     //origin for drawing space and species
  private transient final int[] phaseSize = new int[Space.D];  //array form of width, height
  
 /**
  * Symmetric array of all two-body potentials.  Potentials are associated with species, and each species
  * is assigned a unique index to idenfity it.  Potential2[i][j] is the two-body potential
  * for Species indexed i and j, respectively.  The potential for i=j is merely the one describing the 
  * molecular interactions for that species.
  * 
  * @see Species#speciesIndex
  * @see Potential2
  */
  public Potential2[][] potential2;
  
 /**
  * Array of all one-body potentials.  Potentials are associated with species, and each species
  * is assigned a unique index to idenfity it.  Potential1[i] is the one-body potential
  * for Species indexed i.
  * 
  * @see Species#speciesIndex
  * @see Potential1
  */
  public Potential1[] potential1;
  
 /**
  * Vector of all species
  */
  public transient Vector speciesVector = new Vector(3);
 
  public boolean updatedForces, updatedPotentialEnergy, updatedKineticEnergy;
  public boolean updatedNeighbors, updatedFutureNeighbors;
  
  private double potentialEnergy, kineticEnergy;  // must access with get methods
 
 /**
  * The Space object associated with this Phase.
  */
  Space space;
 
 /**
  * Object used to describe presence and magnitude of constant gravitational acceleration
  */
  public Gravity gravity;
  public boolean noGravity = true;
    
 /**
  * When using periodic boundaries, image molecules near the cell boundaries often have parts that overflow
  * into the central cell.  When the phase is drawn, these "overflow portions" are not normally
  * included in the central image.  Setting this flag to <code>true</code> causes extra drawing
  * to be done so that the overflow portions are properly rendered.  This is particularly helpful
  * to have on when nShells is non-zero.  Default value is <code>false</code>.
  */
  boolean drawOverflowImages = false;
 
 /**
  * First species in the linked list of species in this phase.
  */
   Species firstSpecies;
 
 /**
  * Last species in the linked list of species in this phase.
  */
  Species lastSpecies;
 
 /**
  * First molecule in the linked list of molecules in this phase
  */
//  public Molecule firstMolecule;
  
 /**
  * Last molecule in the linked list of molecules in this phase
  */
//  public Molecule lastMolecule;
 
 /**
  * First atom in the linked list of atoms in this phase
  */
//  public Atom firstAtom;
 
 /**
  * Last atom in the linked list of atoms in this phase
  */
//  public Atom lastAtom;
 
 /**
  * Total number of atoms in this phase
  */
  public int nAtomTotal;
 
 /**
  * Total number of molecules in this phase
  *
  * @see Species#addMolecule
  * @see Species#deleteMolecule
  */
  public int nMoleculeTotal;
 
 /**
  * Image object used for double-buffering
  */
  Image offScreenImage;
 
 /**
  * Graphics object used for double-buffering
  */
  Graphics offScreenGraphics;
  
 /**
  * Initial temperature, in Kelvins
  */
  public double initialTemperature = 300.;
  
  Simulation parentSimulation;
  Meter firstMeter, lastMeter;
  private int nMeters = 0;
  
  private Potential1 p1Null = new P1Null();
  private Potential2 p2IdealGas = new P2IdealGas();
  public Configuration configuration;
    
  public Phase() {
    setLayout(null);
    setSize(300,300);
    setBackground(Color.white);
    updatedForces = updatedPotentialEnergy = updatedKineticEnergy = false;
    updatedNeighbors = updatedFutureNeighbors = false;
    useNeighborList = false;
    nAtomTotal = nMoleculeTotal = 0;
    gravity = new Gravity(0.0);
    noGravity = true;
    add(new ConfigurationSequential());  //default configuration
  }
 
   /********************Eliza's Code**********************/
  /*public void initializeFrame(){  
   
      Button done = new Button("Done");
      done.addMouseListener(new java.awt.event.MouseAdapter()
   	                         {
             	               public void mouseReleased(MouseEvent e)
                    	         {
				                  choice.setVisible(false);
                                  pickConstructor();
                                  Repaint();  
                         	     }
                           	});
    
      choice.setSize(200,200);
	  choice.setBackground(Color.gray);

      done.setBackground(Color.white);
	  done.setBounds(90,150,45,30);
	  
	  configlist.setBackground(Color.white);
	  configlist.addItem("Configuration1");
	  configlist.addItem("Configuration2");
      configlist.addMouseListener(new java.awt.event.MouseAdapter()
                   		        {
                           		  public void mouseReleased(MouseEvent e)
                               	   {
                               	    c = configlist.getSelectedIndex();
                                   }
                               	}); 
    
      mypanel.add(configlist);
	  choice.add(done);
	  choice.add(mypanel);
	  choice.setVisible(true);
  }


  private void pickConstructor(){
    switch(c){
      case 0 : configuration= new Configuration1();
                break;
      case 1 : configuration= new Configuration2();
                break;
     }
  
    C.add(firstSpecies);

  }
  //updates screen after done button is pushed and frame closes
  private void Repaint(){
    repaint();
  }
  */
   /****************End of Eliza's Code*******************/
 
 /**
  * Returns the temperature (in Kelvin) of this phase as computed via the equipartition
  * theorem from the kinetic energy summed over all (atomic) degrees of freedom
  */  
  public double getKineticTemperature() {
    updateKineticEnergy();
    return (2./(double)(nAtomTotal*Space.D))*kineticEnergy*Constants.KE2T;
  }
  
  public void setInitialTemperature(double t) {initialTemperature = t;}
  public double getInitialTemperature() {return initialTemperature;}
  
  public final Atom firstAtom() {
     Molecule m = firstMolecule();
     return (m != null) ? m.firstAtom() : null;
  }
  public final Atom lastAtom() {
    Molecule m = lastMolecule();
    return (m != null) ? m.lastAtom() : null;
  }
  public final Molecule firstMolecule() {
    for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
        Molecule m = s.firstMolecule();
        if(m != null) {return m;}
    }
    return null;
  }
  public final Molecule lastMolecule() {
    for(Species s=lastSpecies; s!=null; s=s.getPreviousSpecies()) {
        Molecule m = s.lastMolecule();
        if(m != null) {return m;}
    }
    return null;
  }
  public final Species firstSpecies() {return firstSpecies;}
  public final Species lastSpecies() {return lastSpecies;}
  
  /** 
   * Resets the class variable Phase.TO_PIXELS so that the x-dimension of this Phase's
   * space exactly fills this Phase's drawing region when space is unscaled
   * setTO_PIXELS = (width of this Phase, in pixels)/(width of Space, in Angstroms)
   *
   * @see #paint
   */
  public void resetTO_PIXELS() {
    if(maySetTO_PIXELS && space!=null) {
        setTO_PIXELS((double)getSize().width/space.getDimensions(0));
    }
  }
  
  /**
   * Sets the value of maySetTO_PIXELS for this phase.  If set to true,
   * resetTO_PIXELS is then called
   *
   * @param b the new value of maySetTO_PIXELS
   * @see #TO_PIXELS
   * @see #resetTO_PIXELS
   */
  public void setMaySetTO_PIXELS(boolean b) {
    maySetTO_PIXELS = b;
    resetTO_PIXELS();
  }
  
 /**
  * @return current value of maySetTO_PIXELS
  * @see #TO_PIXELS
  */
  public boolean getMaySetTO_PIXELS() {return maySetTO_PIXELS;}
 
 /**
  * Sets TO_PIXELS to the given value if maySetTO_PIXELS is <code>true</code>
  */
  public void setTO_PIXELS(double t) {if(maySetTO_PIXELS) {TO_PIXELS = t;}}
 
 /**
  * @return the current value of TO_PIXELS
  */
  public double getTO_PIXELS() {return TO_PIXELS;}
 
 /**
  * Performs a conversion from Angstroms to pixels
  * @param x a length, in simulation units (Angstroms)
  * @return x*TO_PIXELS, the input, converted to pixels
  */
  public static int toPixels(double x) {return (int)(TO_PIXELS*x);}

  public boolean getUseNeighborList() {return useNeighborList;}
  public void setUseNeighborList(boolean b) {useNeighborList = b;}

  public void setDrawBoundingBox(boolean b) {drawBoundingBox = b;}
  public boolean getDrawBoundingBox() {return drawBoundingBox;}
  
  public final double getG() {return gravity.getG();}
  public void setG(double g) {
    gravity.setG(g);
    noGravity = (g == 0.0);
  }
  
 /**
  * @return the current value of imageShells
  */
  public int getImageShells() {return imageShells;}
 
 /**
  * Changes the value of image shells, and calls the setScale method of the Phase's Space
  *
  * @param n the new value of imageShells
  * @see Space#setScale
  */
  public void setImageShells(int n) {
      if(n>=0) {
        imageShells = n;
        if(space != null) {space.setScale(nominalScale,imageShells);}
      }
  }
  
  // STOPPED commenting here
  
  
  public double getNominalScale() {return nominalScale;}
  public void setNominalScale(double s) {
      if(s>0) {
        nominalScale = s;
        if(space != null) {space.setScale(nominalScale,imageShells);}
      }
  }
    
  //Override superclass methods for changing size so that TO_PIXELS is reset with any size change  
  // this setBound is ultimately called by all other setSize, setBounds methods
  public void setBounds(int x, int y, int width, int height) {
    super.setBounds(x,y,width,height);
    phaseSize[0] = width;
    phaseSize[1] = height;
    resetTO_PIXELS();
    if(space != null) {space.resetOrigins(imageShells);}
  }
  
  public int[] getPhaseSize() {return phaseSize;}
   
  // end of size-change method overrides

    public void add(Species species) {
        super.add(species);
        species.parentPhase = this;
        species.configurationMolecule.initializeCoordinates();
        if(space != null) species.initializeSpecies(this);
        configuration.add(species);
        speciesVector.addElement(species);
        if(lastSpecies != null) {lastSpecies.setNextSpecies(species);}
        else {firstSpecies = species;}
        lastSpecies = species;
        for(Molecule m=species.firstMolecule(); m!=null; m=m.getNextMolecule()) {nMoleculeTotal++;}
        for(Atom a=species.firstAtom(); a!=null; a=a.getNextAtom()) {nAtomTotal++;}
        if(species.getSpeciesIndex() > speciesCount-1) {setSpeciesCount(species.getSpeciesIndex()+1);}
    }
    
    /* Resizes potential arrays, keeping all elements already filled in, and
       setting to p1Null or p2IdealGas the newly formed elements
       */
    private void setSpeciesCount(int n) {
        Potential1 p1[] = new Potential1[n];
        Potential2 p2[][] = new Potential2[n][n];
        for(int i=0; i<speciesCount; i++) {
            p1[i] = potential1[i];
            p2[i][i] = potential2[i][i];
            for(int j=i+1; j<speciesCount; j++) {
                p2[i][j] = p2[j][i] = potential2[i][j];
                
            }
        }
        for(int i=speciesCount; i<n; i++) {
            p1[i] = p1Null;
            p2[i][i] = p2IdealGas;
            for(int j=0; j<n; j++) {
                p2[i][j] = p2[j][i] = p2IdealGas;
            }
        }
        potential1 = p1;
        potential2 = p2;
        speciesCount = n;
    }
    
/*    public final void setFirstMolecule() {
        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
            if(s.firstMolecule != null) {
//                firstMolecule = s.firstMolecule;
//                firstAtom = s.firstAtom();
                return;
            }
        }
    }
    
    
    public final void setLastMolecule() {
        for(Species s=lastSpecies; s!=null; s=s.getPreviousSpecies()) {
            if(s.lastMolecule != null) {
//                lastMolecule = s.lastMolecule;
//                lastAtom = s.lastAtom();
                return;
            }
        }
    }
*/    
  
     public void add(Space space) {
        super.add(space);
        this.space = space;
        space.setParentPhase(this);
        int w = getSize().width;
        int h = getSize().height;
        space.setBounds(0,0,w,h);
        space.setScale(nominalScale, imageShells);
        resetTO_PIXELS();
        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {s.initializeSpecies(this);}
        for(int i = 0; i < speciesCount; i++) {
            if(potential1[i] != null) potential1[i].setSpace(space);
            for(int j = i; j < speciesCount; j++) {
                if(potential2[i][j] != null) potential2[i][j].setSpace(space);
            }
        }   
        // Instantiate some things for double-buffering here because cannot be done in Phase's constructor
        if(offScreenImage==null) {
            offScreenImage = this.createImage(getSize().width,getSize().height);
            offScreenGraphics = offScreenImage.getGraphics();
        }
    }
    
    public void add(Potential1 p1) {
        super.add(p1);
        p1.setPhase(this);
        if(space != null) p1.setSpace(space);
        if(p1.speciesIndex+1 > speciesCount) {setSpeciesCount(p1.speciesIndex+1);}
        this.potential1[p1.speciesIndex] = p1;
    }
    
    public void add(Potential2 p2) {
        super.add(p2);
        p2.setPhase(this);
        if(space != null) p2.setSpace(space);
        int idx = Math.max(p2.species1Index,p2.species2Index);
        if(idx+1 > speciesCount) {setSpeciesCount(idx+1);}
        this.potential2[p2.species1Index][p2.species2Index] = p2;
        this.potential2[p2.species2Index][p2.species1Index] = p2;
/*        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
            if(p2.species1Index == s.getSpeciesIndex() || p2.species2Index == s.getSpeciesIndex()) {
                double d = 0.25*p2.skinThickness*p2.skinThickness;
                if(d < s.getNeighborUpdateSquareDisplacement()) {
                    s.setNeighborUpdateSquareDisplacement(d);
                }
            }
          }
        */
    }
    
    public void add(Configuration c){
        c.parentPhase = this;
        configuration = c;
        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
            configuration.add(s);
        }
    }
    
	public void add(Meter m) {
	    if(firstMeter == null) {firstMeter = m;}
	    if(lastMeter != null) {lastMeter.setNextMeter(m);}
	    lastMeter = m;
	    nMeters++;
	    m.phase = this;
	    m.initialize();
	    if(parentSimulation.haveIntegrator()) {
	        parentSimulation.controller.integrator.addIntegrationIntervalListener(m);
	    }
	}
	
	// Returns ith meter in linked list of meters, with i=0 being the first meter
	public Meter getMeter(int i) {
	    if(i >= nMeters) {return null;}
	    Meter m = firstMeter;
        for(int j=i; --j>=0; ) {m = m.getNextMeter();}  //get ith meter in list
        return m;
    }

    // Updates neighbor list for all molecules
    // Only preceding molecules are included in a given molecule's neighbor list
    public void updateNeighbors() {
        for(Molecule m1=firstMolecule(); m1!=null; m1=m1.getNextMolecule()) {
            m1.clearNeighborList();
            for(Molecule m2=firstMolecule(); m2!=m1; m2=m2.getNextMolecule()) {
                if(potential2[m1.getSpeciesIndex()][m2.getSpeciesIndex()].isNeighbor(m1,m2)) {
                    m1.addNeighbor(m2);
                }
            }
        }
        updatedNeighbors = true;
    }
        
    public void updateForces() {
/*      if(updatedForces) {return;}
      potentialEnergy = 0.0;
      for(Molecule m1=firstElement; m1!=null; m1=m1.getNextMolecule()) {
        m1.zeroForce();
        Molecule m2;
//        for(Molecule m2=firstElement; m2!=m1; m2=m2.getNext()) {
        for(Enumeration enum=m1.getNeighborList(); enum.hasMoreElements();) {
           m2 = (Molecule)enum.nextElement();
           PairInteraction pair = potential[m1.getSpeciesIndex()][m2.getSpeciesIndex()].computePairInteraction(m1,m2);
           m2.addForce(pair.force);
           m1.subtractForce(pair.force);
           potentialEnergy += pair.energy;
        }
      }
      updatedForces = true;
      updatedPotentialEnergy = true;
*/    }
    
    public void updateKineticEnergy() {
      kineticEnergy = 0.0;
      for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
        kineticEnergy += s.kineticEnergy();
      }
      updatedKineticEnergy = true;
    }
 
 // Works only for one species
    public boolean overlap(Atom atom, double energy) {
      energy = 0.0;
      for(Atom a=firstAtom(); a!=atom; a=a.getNextAtom()) {
         double u = 0.0;
         Potential pot = potential2[0][0].getPotential(atom,a);
         if(pot.overlap(atom,a,u)) {return true;}
         else {energy += u;}
      }
      for(Atom a=atom.getNextAtom(); a!=null; a=a.getNextAtom()) {
         double u = 0.0;
         Potential pot = potential2[0][0].getPotential(atom,a);
         if(pot.overlap(atom,a,u)) {return true;}
         else {energy += u;}
      }      
      return false;
    }
        
    public void updatePotentialEnergy() {
/*      potentialEnergy = 0.0;
      for(Molecule m1=firstMolecule; m1!=null; m1=m1.getNextMolecule()) {
        for(Molecule m2=firstMolecule; m2!=m1; m2=m2.getNextMolecule()) {
          PairInteraction pair = potential[m1.getSpeciesIndex()][m2.getSpeciesIndex()].computePairInteraction(m1,m2);
          potentialEnergy += pair.energy;
        }
      }
      updatedPotentialEnergy = true;
*/    }

//needs work for efficiency and multiatomics   
    public double potentialEnergy() {  
        double energy = 0.0;
        for(Molecule m1=firstMolecule(); m1!=null; m1=m1.getNextMolecule()) {
           energy += m1.potentialEnergy();
        }
        return (0.5*energy);
      }
        

    public double getKineticEnergy() {
        if(!updatedKineticEnergy) {updateKineticEnergy();}
        return kineticEnergy;
    }
    
    public double getPotentialEnergy() {
        if(!updatedPotentialEnergy) {updatePotentialEnergy();}
        return potentialEnergy;
    }
    
    public double getTotalEnergy() {
        return getKineticEnergy() + getPotentialEnergy();
    }
      
  public final boolean getDrawOverflowImages() {return drawOverflowImages;}
  public final void setDrawOverflowImages(boolean b) {drawOverflowImages = b;}

 /** 
  * paint is the method that handles the drawing of the phase to the screen.
  * Several variables and conditions affect how the image is drawn.  First,
  * the class variable <code>TO_PIXELS</code> performs the conversion between simulation
  * length units (Angstroms) and pixels.  The default value is 300 pixels/Angstrom
  * reflecting the default size of the phase (300 pixels by 300 pixels) and the
  * default length scale (selected so that the simulation volume is of unit 
  * length).  The size of the phase (in simulation units) is held in <code>space.dimensions[]</code>.
  * Two quantities can further affect the size of the drawn image of the phase. 
  * The variable <code>nominalScale</code> is a multiplicative factor that directly
  * scales up or down the size of the image; scaling of the image is also performed
  * whenever shells of periodic images are drawn.  Scaling is performed automatically
  * to permit the central image and all of the specified periodic images to fit in
  * the drawing of the phase.  The number of image shells, together with the nominalScale,
  * are taken by <code>space</code> to determine the overall scaling of the drawn image.
  *
  * Painting is done with double buffering.  First a solid rectangle of the background
  * color is drawn to an off-screen graphic.  Then the origin of the drawn image is
  * determined from the size of the drawn image and the size of the phase:
  * origin[i] = 0.5*(phaseSize[i] - spaceSize[i]).  A gray line is drawn around
  * the phase boundary if <code>drawBoundingBox</code> is <code>true</code>.
  * A loop is then taken over all species, first passing the species to space.repositionMolecules
  * to enforce (for example) periodic boundaries, then invoking the draw method
  * of the species.  The draw method of <code>space</code> is then invoked.
  * If imageShells is non-zero, the getImageOrigins method of space is invoked to
  * determine the origins for all images, and a replica of the just-completed image
  * of the central cell is copied to all of the periodic images.  The complete
  * off-screen image is then transferred to the graphics object g.  
  *
  * Note that handling
  * of overflowImages (parts of neighboring periodic images that spill into the
  * central image, and which must be rendered separately) is performed by the 
  * species draw method.
  *
  * @param g The graphic object to which the image of the phase is drawn
  * @see Space
  * @see Species
  */
  public void paint(Graphics g) {
    if(Beans.isDesignTime()){
        g.setColor(getBackground());
        g.drawRect(0,0,getSize().width-1,getSize().height-1);
        g.drawRect(1,1,getSize().width-3,getSize().height-3);
        paintComponents(g);
    }
    else {
        int w = getSize().width;
        int h = getSize().height;
        offScreenGraphics.setColor(getBackground());
        offScreenGraphics.fillRect(0,0,w,h);
        if(drawBoundingBox) {
            offScreenGraphics.setColor(Color.gray);
            offScreenGraphics.drawRect(0,0,w-1,h-1);
            }
        origin = space.getCentralOrigin();
        for(int i = 0; i < speciesVector.size(); i++) {
            Species species = (Species)speciesVector.elementAt(i);
            if(species.firstAtom() == null) {continue;}
            space.repositionMolecules(species);
            species.draw(offScreenGraphics, origin, space.getScale());
            }
        space.draw(offScreenGraphics, origin);
        origin = space.getCopyOrigin();
        if(imageShells > 0) {
            int[][] origins = space.getImageOrigins(imageShells);
            int[] spaceSize = space.getDrawSize();
            for(int i=0; i<origins.length; i++) {
            
    /*        // brute-force way to make periodic images (doesn't work well)
                for(int j = 0; j < speciesVector.size(); j++) {
                    Species species = (Species)speciesVector.elementAt(j);
//                    space.repositionElements(species);
                    species.draw(offScreenGraphics, origins[i], space.getScale());
                }
                space.draw(offScreenGraphics, origins[i]);
             // */   
             // elegant way
                offScreenGraphics.copyArea(origin[0],origin[1],spaceSize[0],spaceSize[1],origins[i][0],origins[i][1]);
            }
        }
        g.drawImage(offScreenImage,0,0,this);
    }
  }    
}
    