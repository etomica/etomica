package simulate;
import java.awt.*;//Because this class subclasses Component.
import java.io.*;//for Serializable and Stream.
import java.util.*; //for type Vector
import java.beans.*;//for PropertyChangeSupport and other beans convenience classes.
import java.util.*; //for neighborList

public final class Phase extends Container {

  public static double TO_PIXELS = 300.;
  
  public  boolean useNeighborList;
  protected transient Vector listenerspecies=null;
  private Vector phaseIntegratorListeners = new Vector(3);
  private boolean drawBoundingBox = false;
  private boolean maySetTO_PIXELS = true;
  private int imageShells = 0;
  private double nominalScale = 1.0;
  int speciesCount=0;
  private transient final int[] origin = new int[Space.D];     //origin for drawing space and species
  private transient final int[] phaseSize = new int[Space.D];  //array form of width, height
  public Potential2[][] potential2;
  public Potential1[] potential1;
  public transient Vector speciesVector = new Vector(3);
  public boolean updatedForces, updatedPotentialEnergy, updatedKineticEnergy;
  public boolean updatedNeighbors, updatedFutureNeighbors;
  private double potentialEnergy, kineticEnergy;  // must access with get methods
  Space space;
  boolean drawOverflowImages = false;
  public Species firstSpecies, lastSpecies;
  public Molecule firstMolecule, lastMolecule;
  public Atom firstAtom, lastAtom;
  public int nAtomTotal, nMoleculeTotal;  //totals in this phase
  Image offScreenImage;
  public double initialTemperature = 300.;
  Graphics offScreenGraphics;
  
  public Phase() {
    setLayout(null);
    setSize(300,300);
    setBackground(Color.gray);
    updatedForces = updatedPotentialEnergy = updatedKineticEnergy = false;
    updatedNeighbors = updatedFutureNeighbors = false;
    useNeighborList = false;
    nAtomTotal = nMoleculeTotal = 0;
  }
  
  public double getKineticTemperature() {
    updateKineticEnergy();
    return (2./(double)(nAtomTotal*Space.D))*kineticEnergy*Constants.KE2T;
  }
  public void setInitialTemperature(double t) {initialTemperature = t;}
  public double getInitialTemperature() {return initialTemperature;}
  
  /** Resets the class variable Phase.TO_PIXELS so that the Phase's
      space exactly fills the Phase's drawing region when space is unscaled
      */
  public void resetTO_PIXELS() {
    if(maySetTO_PIXELS && space!=null) {
        setTO_PIXELS((double)getSize().width/space.getDimensions(0));
    }
  }
  
  public void setMaySetTO_PIXELS(boolean b) {
    maySetTO_PIXELS = b;
    resetTO_PIXELS();
  }
  public boolean getMaySetTO_PIXELS() {return maySetTO_PIXELS;}
  
  public void setTO_PIXELS(double t) {TO_PIXELS = t;}
  public double getTO_PIXELS() {return TO_PIXELS;}
        
  public static int toPixels(double x) {return (int)(TO_PIXELS*x);}

  public boolean getUseNeighborList() {return useNeighborList;}
  public void setUseNeighborList(boolean b) {useNeighborList = b;}

  /** set and get flag that determines whether a box is drawn around 
      Phase's drawing region
      */
  public void setDrawBoundingBox(boolean b) {drawBoundingBox = b;}
  public boolean getDrawBoundingBox() {return drawBoundingBox;}
  
  public int getImageShells() {return imageShells;}
  public void setImageShells(int n) {
      if(n>=0) {
        imageShells = n;
        if(space != null) {space.setScale(nominalScale,imageShells);}
      }
  }
  
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
  }
  
  public int[] getPhaseSize() {return phaseSize;}
   
  // end of size-change method overrides
 
  public synchronized void addPhaseIntegratorListener(PhaseIntegratorListener pil) {
    phaseIntegratorListeners.addElement(pil);
    firePhaseIntegratorEvent(new PhaseIntegratorEvent(this));
  }

  public synchronized void removePhaseIntegratorListener(PhaseIntegratorListener pil) {
    phaseIntegratorListeners.removeElement(pil);
  }

  public void firePhaseIntegratorEvent(PhaseIntegratorEvent pie) {
    Vector currentListeners = null;
    synchronized(this){
        currentListeners = (Vector)phaseIntegratorListeners.clone();
    }
    for(int i = 0; i < currentListeners.size(); i++) {
        PhaseIntegratorListener listener = (PhaseIntegratorListener)currentListeners.elementAt(i);
        listener.phaseIntegratorNotify(pie);
    }
  }
 
    public void add(Species species) {
        super.add(species);
        species.initializeSpecies(this);
        speciesVector.addElement(species);
        if(speciesCount > 0) {lastSpecies.setNextSpecies(species);}
        else {firstSpecies = species;}
        lastSpecies = species;
        setFirstMolecule();
        setLastMolecule();
        speciesCount++;
        for(Molecule m=species.firstMolecule; m!=null; m=m.getNextMolecule()) {nMoleculeTotal++;}
        for(Atom a=species.firstAtom; a!=null; a=a.getNextAtom()) {nAtomTotal++;}
        potential1 = new Potential1[speciesCount];
        potential2 = new Potential2[speciesCount][speciesCount];
        potential1[0] = new P1Null();
        potential2[0][0] = new P2IdealGas();
        for(int i = 0; i < speciesCount; i++) {
            potential1[i] = potential1[0];
            potential2[i][i] = potential2[0][0];
            for(int j = i + 1; j < speciesCount; j++) {
                potential2[i][j] = potential2[0][0];
                potential2[j][i] = potential2[i][j];
            }
        }   
    }
    
    public final void setFirstMolecule() {
        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
            if(s.firstMolecule != null) {
                firstMolecule = s.firstMolecule;
                firstAtom = s.firstAtom;
                return;
            }
        }
    }
    
    public final void setLastMolecule() {
        for(Species s=lastSpecies; s!=null; s=s.getPreviousSpecies()) {
            if(s.lastMolecule != null) {
                lastMolecule = s.lastMolecule;
                lastAtom = s.lastAtom;
            }
        }
    }
  
     public void add(Space space) {
        super.add(space);
        this.space = space;
        space.setParentPhase(this);
        int w = getSize().width;
        int h = getSize().height;
        space.setBounds(0,0,w,h);
        space.setScale(nominalScale, imageShells);
        resetTO_PIXELS();
        
        // Instantiate some things for double-buffering here because cannot be done in Phase's constructor
        if(offScreenImage==null) {
            offScreenImage = this.createImage(getSize().width,getSize().height);
            offScreenGraphics = offScreenImage.getGraphics();
        }
    }
    
    public void add(Potential1 p1) {
        super.add(p1);
        p1.setSpace(this.space);
        this.potential1[p1.speciesIndex] = p1;
    }
    
    public void add(Potential2 p2) {
        super.add(p2);
        p2.setSpace(this.space);
        this.potential2[p2.species1Index][p2.species2Index] = p2;
        this.potential2[p2.species2Index][p2.species1Index] = p2;
        for(Species s=firstSpecies; s!=null; s=s.getNextSpecies()) {
            if(p2.species1Index == s.getSpeciesIndex() || p2.species2Index == s.getSpeciesIndex()) {
                double d = 0.25*p2.skinThickness*p2.skinThickness;
                if(d < s.getNeighborUpdateSquareDisplacement()) {
                    s.setNeighborUpdateSquareDisplacement(d);
                }
            }
        }
    }
    
    // Updates neighbor list for all molecules
    // Only preceding molecules are included in a given molecule's neighbor list
    public void updateNeighbors() {
        for(Molecule m1=firstMolecule; m1!=null; m1=m1.getNextMolecule()) {
            m1.clearNeighborList();
            for(Molecule m2=firstMolecule; m2!=m1; m2=m2.getNextMolecule()) {
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
      for(Atom a=firstAtom; a!=atom; a=a.getNextAtom()) {
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
        for(Molecule m1=firstMolecule; m1!=null; m1=m1.getNextMolecule()) {
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
        offScreenGraphics.setColor(Color.white);
        offScreenGraphics.fillRect(0,0,w,h);
        int[] spaceSize = space.getDrawSize();
        Space.uEa1T_v1Mv2_(origin,0.5,phaseSize,spaceSize);
        if(drawBoundingBox) {
            offScreenGraphics.setColor(Color.gray);
            offScreenGraphics.drawRect(0,0,w-1,h-1);
            }
        for(int i = 0; i < speciesVector.size(); i++) {
            Species species = (Species)speciesVector.elementAt(i);
            if(species.firstAtom == null) {continue;}
            space.repositionMolecules(species);
            species.draw(offScreenGraphics, origin, space.getScale());
            }
        space.draw(offScreenGraphics, origin);
        if(imageShells > 0) {
            int[][] origins = space.getImageOrigins(imageShells);
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
  
  
/* Don't know where this came from
	class SymContainer extends java.awt.event.ContainerAdapter
	{
		public void componentRemoved(java.awt.event.ContainerEvent event)
		{
			Object object = event.getSource();
			if (object == Phase.this)
				Phase_ComponentRemoved(event);
		}
	}

	void Phase_ComponentRemoved(java.awt.event.ContainerEvent event)
	{
		// to do: code goes here.
	}
	*/
}
    