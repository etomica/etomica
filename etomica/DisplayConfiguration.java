package simulate;
import java.awt.*;
import java.beans.Beans;

    public class DisplayConfiguration extends Display {

	Simulation parentSimulation;
    int pixels = 200;
    Image offScreen;
    Graphics osg;
    int updateInterval;
    int iieCount;
    Component displayTool = null;

 /**
  * Flag specifying whether a line tracing the boundary of the phase should be 
  * included when drawing the phase to the screen.
  * Default value is <code>false</code>
  */
  private boolean drawBoundingBox = false;
  
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
  
  private transient int[] origin;     //origin for drawing space and species
  private transient final int[] phaseSize = new int[Space.D];  //array form of width, height
 /**
  * When using periodic boundaries, image molecules near the cell boundaries often have parts that overflow
  * into the central cell.  When the phase is drawn, these "overflow portions" are not normally
  * included in the central image.  Setting this flag to <code>true</code> causes extra drawing
  * to be done so that the overflow portions are properly rendered.  This is particularly helpful
  * to have on when nShells is non-zero.  Default value is <code>false</code>.
  */
  boolean drawOverflowImages = false;
 

    public DisplayConfiguration () {
        setSize(pixels, pixels);
        setBackground(Color.white);
	    setUpdateInterval(1);
    }

    public void paint(Graphics g) {
      if(Beans.isDesignTime()) {
        g.setColor(Color.red);
        g.drawRect(0,0,getSize().width-1,getSize().height-1);
        g.drawRect(1,1,getSize().width-3,getSize().height-3);
      } 
      createOffScreen();
      doPaint(osg);
      g.drawImage(offScreen, 0, 0, null);
    }

  public final boolean getDrawOverflowImages() {return drawOverflowImages;}
  public final void setDrawOverflowImages(boolean b) {drawOverflowImages = b;}

  public double getNominalScale() {return nominalScale;}
  public void setNominalScale(double s) {
      if(s>0) {
        nominalScale = s;
        if(phase.space != null) {phase.space.setScale(nominalScale,imageShells);}
      }
  }
    
  //Override superclass methods for changing size so that TO_PIXELS is reset with any size change  
  // this setBound is ultimately called by all other setSize, setBounds methods
  public void setBounds(int x, int y, int width, int height) {
    super.setBounds(x,y,width,height);
    phaseSize[0] = width;
    phaseSize[1] = height;
    phase.resetTO_PIXELS();
    if(phase.space != null) {phase.space.resetOrigins(imageShells);}
  }
  
  public int[] getPhaseSize() {return phaseSize;}
   
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
        if(phase.space != null) {phase.space.setScale(nominalScale,imageShells);}
      }
  }
  
  public void setDrawBoundingBox(boolean b) {drawBoundingBox = b;}
  public boolean getDrawBoundingBox() {return drawBoundingBox;}

  public void doUpdate() {;}
    
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
  public void doPaint(Graphics g) {
    if(Beans.isDesignTime()){
        g.setColor(getBackground());
        g.drawRect(0,0,getSize().width-1,getSize().height-1);
        g.drawRect(1,1,getSize().width-3,getSize().height-3);
    }
    else {
        if(phase == null) {return;}
        int w = getSize().width;
        int h = getSize().height;
        g.setColor(getBackground());
        g.fillRect(0,0,w,h);
        if(drawBoundingBox) {
            g.setColor(Color.gray);
            g.drawRect(0,0,w-1,h-1);
            }
        Space space = phase.space;
        if(space == null) {return;}
        origin = space.getCentralOrigin();
        for(Species s=phase.firstSpecies(); s!=null; s=s.getNextSpecies()) {
            if(s.firstAtom() == null) {continue;}
            space.repositionMolecules(s);
            s.draw(g, origin, space.getScale());
            }
        space.draw(g, origin);
        origin = space.getCopyOrigin();
        if(imageShells > 0) {
            int[][] origins = space.getImageOrigins(imageShells);
            int[] spaceSize = space.getDrawSize();
            for(int i=0; i<origins.length; i++) {
                g.copyArea(origin[0],origin[1],spaceSize[0],spaceSize[1],origins[i][0],origins[i][1]);
            }
        }
    }
  }    
    
}
