package simulate;
import java.awt.*;
import java.beans.Beans;

    public class DisplayConfiguration extends Display {

    public static final int LEFT = -1;   //Class variables to code for alignment of drawn image within display region
    public static final int CENTER = 0;
    public static final int RIGHT = +1;
    public static final int TOP = -1;
    public static final int BOTTOM = +1;

    public final int align[] = new int[Space.D];
    
 /**
  * Size of drawing region of central image, in pixels
  *
  * @see #computeDrawSize
  */
    protected final int[] drawSize = new int[Space.D];
  
 /**
  * Flag specifying whether a line tracing the boundary of the display should be drawn
  * Default value is <code>true</code>
  */
  private boolean drawBoundingBox = true;
  
 /**
  * Flag specifying whether a line tracing the boundary of the space should be drawn
  * Default value is <code>false</code>
  */
  private boolean drawSpaceFrame = false;
  
 /**
  * Number of periodic-image shells to be drawn when drawing this phase to the
  * screen.  Default value is 0.
  *
  * @see #paint
  */
  private int imageShells = 0;
 
 /**
  * Factor used to scale the size of the image. May be used
  * to scale up or down the image within one phase without affecting those
  * in other displays.  Default value is 1.0.
  *
  * @see Phase#paint
  */
    protected double scale = 1.0;
      
   /**
    * Coordinate origin for central image
    */
    protected final int[] centralOrigin = new int[Space.D];

 /**
  * When using periodic boundaries, image molecules near the cell boundaries often have parts that overflow
  * into the central cell.  When the phase is drawn, these "overflow portions" are not normally
  * included in the central image.  Setting this flag to <code>true</code> causes extra drawing
  * to be done so that the overflow portions are properly rendered.  This is particularly helpful
  * to have on when imageShells is non-zero.  Default value is <code>false</code>.
  */
  public static boolean DRAW_OVERFLOW = false;
 

    public DisplayConfiguration () {
        super();
	    Space.uEa1(align,CENTER);   //default alignment is centered in x and y directions
    }
    
    public void setAlign(int i, int value) {
        align[i] = value;
    }
    public int getAlign(int i) {return align[i];}

  public final boolean getDRAW_OVERFLOW() {return DRAW_OVERFLOW;}
  public final void setDRAW_OVERFLOW(boolean b) {DRAW_OVERFLOW = b;}

  public double getScale() {return scale;}
  public void setScale(double s) {
      if(s>0) {
        scale = s;
      }
  }
    
  //Override superclass methods for changing size so that TO_PIXELS is reset with any size change  
  // this setBound is ultimately called by all other setSize, setBounds methods
  public void setBounds(int x, int y, int width, int height) {
    if(getBounds().width * getBounds().height != 0) {  //reset scale based on larger size change
        double ratio1 = (double)width/(double)getBounds().width;
        double ratio2 = (double)height/(double)getBounds().height;
        double factor = Math.min(ratio1, ratio2);
//        double factor = (Math.abs(Math.log(ratio1)) > Math.abs(Math.log(ratio2))) ? ratio1 : ratio2;
        scale *= factor;
    }
    super.setBounds(x,y,width,height);
//    if(phase != null) phase.resetTO_PIXELS();
  }
   
 /**
  * @return the current value of imageShells
  */
  public int getImageShells() {return imageShells;}
 
 /**
  * Changes the value of image shells, and increases/decreases scale accordingly
  *
  * @param n the new value of imageShells
  */
  public void setImageShells(int n) {
      if(n>=0) {
        scale *= (double)(2*imageShells+1)/(double)(2*n+1);
        imageShells = n;
      }
  }
  
    protected int computeOrigin(int align, int drawSize, int size) {
        switch(align) {
            case   LEFT: return 0;    //same as TOP
            case CENTER: return (size-drawSize)/2;
            case  RIGHT: return size-drawSize; //same as BOTTOM
            default: return 0;
        }
    }
        
  public void setDrawBoundingBox(boolean b) {drawBoundingBox = b;}
  public boolean getDrawBoundingBox() {return drawBoundingBox;}
  public void setDrawSpaceFrame(boolean b) {drawSpaceFrame = b;}
  public boolean getDrawSpaceFrame() {return drawSpaceFrame;}

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
  public void doPaint(Graphics g) {  //specific to 2-D
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
        double toPixels = scale*Phase.TO_PIXELS;
        Space.uEa1Tv1(drawSize,toPixels,space.dimensions);
        centralOrigin[0] = computeOrigin(align[0],drawSize[0],w);
        centralOrigin[1] = computeOrigin(align[1],drawSize[1],h);
        for(Species s=phase.firstSpecies(); s!=null; s=s.getNextSpecies()) {
            if(s.firstAtom() == null) {continue;}
            s.draw(g, centralOrigin, scale);
        }
        if(drawSpaceFrame) {space.drawFrame(g, centralOrigin, scale);}
        if(imageShells > 0) {
            double[][] origins = space.imageOrigins(imageShells);  //more efficient to save rather than recompute each time
            for(int i=0; i<origins.length; i++) {
                g.copyArea(centralOrigin[0],centralOrigin[1],drawSize[0],drawSize[1],(int)(toPixels*origins[i][0]),(int)(toPixels*origins[i][1]));
            }
        }
  }        
}
