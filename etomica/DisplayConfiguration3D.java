package simulate;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.beans.Beans;

import org.openscience.jmol.*;
import javax.swing.JSlider;
import javax.swing.*;
import javax.swing.event.*;

/**
 * Displays configuration of molecules
 * Assumes 2-dimensional system
 */
public class DisplayConfiguration3D extends Display {
        
 // ChemFrame cf;  
 // JSlider xs,ys,zs;
 // static Jmol myJmolPanel = new Jmol();

    public static double SIM2PIXELS = 300.;  //Conversion factor from simulation units to display pixels    

    public static final int LEFT = -1;   //Class variables to code for alignment of drawn image within display region
    public static final int CENTER = 0;
    public static final int RIGHT = +1;
    public static final int TOP = -1;
    public static final int BOTTOM = +1;
    private final int D = 2;
    Atom.Iterator downList;
    Atom.Iterator upList;

    public final int[] align = new int[D];
    
 /**
  * Size of drawing region of central image, in pixels
  *
  * @see #computeDrawSize
  */
    protected final int[] drawSize = new int[D];
  
 /**
  * Flag specifying whether a line tracing the boundary of the display should be drawn
  * Default value is <code>true</code>
  */
  private boolean drawBoundingBox = true;
  
 /**
  * Flag specifying whether a line tracing the boundary of the space should be drawn
  * Default value is <code>false</code>
  */
  private boolean drawPhase = false;
  
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
    protected final int[] centralOrigin = new int[D];

 /**
  * When using periodic boundaries, image molecules near the cell boundaries often have parts that overflow
  * into the central cell.  When the phase is drawn, these "overflow portions" are not normally
  * included in the central image.  Setting this flag to <code>true</code> causes extra drawing
  * to be done so that the overflow portions are properly rendered.  This is particularly helpful
  * to have on when imageShells is non-zero.  Default value is <code>false</code>.
  */
  public static boolean DRAW_OVERFLOW = false;
  
  private boolean writeScale = false;  //flag to indicate if value of scale should be superimposed on image
  private TextField scaleText = new TextField();
  private Font font = new Font("sansserif", Font.PLAIN, 10);
//  private int annotationHeight = font.getFontMetrics().getHeight();
  private int annotationHeight = 12;
  private Phase phase2D;
  JmolDisplay display;

    public DisplayConfiguration3D () {
        super();
        align[0] = align[1] = CENTER;
        scaleText.setVisible(true);
        scaleText.setEditable(false);
        scaleText.setBounds(0,0,100,50);
//        JmolDisplay display = new JmolDisplay(this); 
    }
    
    public void setAlign(int i, int value) {
        align[i] = value;
    }
    public int getAlign(int i) {return align[i];}

  public void setPhase(Phase p) {  //2D needed to manipulate dimensions array directly
    super.setPhase(p); phase2D = p;
    downList = p.iterator.makeAtomIteratorDown();
    upList = p.iterator.makeAtomIteratorUp();
    downList.reset(p.firstAtom());
    upList.reset(p.firstAtom());
    display = new JmolDisplay(this); 
    }  
 
    protected int computeOrigin(int align, int drawSize, int size) {
        switch(align) {
            case   LEFT: return 0;    //same as TOP
            case CENTER: return (size-drawSize)/2;
            case  RIGHT: return size-drawSize; //same as BOTTOM
            default: return 0;
        }
    }
        
  public void setDrawPhase(boolean b) {drawPhase = b;}
  public boolean getDrawPhase() {return drawPhase;}

  public void doUpdate() {;}
  
  public void doPaint(Graphics g) {  //specific to 2-D
        if(!isVisible() || phase2D == null) {return;}
 /*       int w = getSize().width;
        int h = getSize().height;
        g.setColor(getBackground());
        g.fillRect(0,0,w,h);
        if(drawBoundingBox) {
            g.setColor(Color.gray);
            g.drawRect(0,0,w-1,h-1);
            }
        double toPixels = scale*SIM2PIXELS;
        drawSize[0] = (int)(toPixels*phase2D.boundary().dimensions().component(0));
        drawSize[1] = (int)(toPixels*phase2D.boundary().dimensions().component(1));
        centralOrigin[0] = computeOrigin(align[0],drawSize[0],w);
        centralOrigin[1] = computeOrigin(align[1],drawSize[1],h);
        for(Species.Agent s=phase2D.firstSpecies(); s!=null; s=s.nextSpecies()) {
            if(s.firstAtom() == null) {continue;}
            s.draw(g, centralOrigin, scale);
        }
        if(drawPhase) {phase2D.paint(g, centralOrigin, scale);}
        if(imageShells > 0) {
            double[][] origins = phase2D.boundary().imageOrigins(imageShells);  //more efficient to save rather than recompute each time
            for(int i=0; i<origins.length; i++) {
                g.copyArea(centralOrigin[0],centralOrigin[1],drawSize[0],drawSize[1],(int)(toPixels*origins[i][0]),(int)(toPixels*origins[i][1]));
            }
        }
        if(writeScale) {
            g.setColor(Color.lightGray);
            g.fillRect(0,getSize().height-annotationHeight,getSize().width,annotationHeight);
            g.setColor(Color.black);
            g.setFont(font);
            g.drawString("Scale: "+Integer.toString((int)(100*scale))+"%", 0, getSize().height-3);
        }
 */
        int atomCount = 0;
/*        upList.reset();
        downList.reset();
        Atom dummy = downList.next();
        while(upList.hasNext()) {
            Atom a = upList.next();
            display.cf.vert[atomCount*3] = (float)a.coordinate().position(0);
            display.cf.vert[atomCount*3+1] = (float)a.coordinate().position(1);
            atomCount++;
        }
        while(downList.hasNext()) {
            Atom a = downList.next();
            display.cf.vert[atomCount*3] = (float)a.coordinate().position(0);
            display.cf.vert[atomCount*3+1] = (float)a.coordinate().position(1);
            atomCount++;
        }
*/            
        for(Atom a=phase2D.firstAtom(); a!=null; a=a.nextAtom()) {
            phase2D.boundary().centralImage(a.coordinate.position());        //move atom to central image
            display.cf.vert[atomCount*3] = (float)a.coordinate().position(0);
            display.cf.vert[atomCount*3+1] = (float)a.coordinate().position(1);
            display.cf.vert[atomCount*3+2] = (float)a.coordinate().position(2);
 //           System.out.println( display.cf.vert[atomCount*3]);
 //           System.out.println( display.cf.vert[atomCount*3+1]);
            atomCount++;
            
        }
//        System.out.println(atomCount);
            display.myJmolPanel.display.repaint();   
  }
  
  public static class JmolDisplay extends JFrame
    {
    ChemFrame cf;  
    JSlider xs,ys,zs;
    static Jmol myJmolPanel = new Jmol();

    public JmolDisplay(DisplayConfiguration3D d3D){
        super ("Jmoldisplay");
        cf.setShowBonds(false);
           
        Jmol.beguine1( (JFrame) this);
        /// * use Jmol.beguine1() instead of myJmolPanel.beguine1()
        /// 'this' refers to JmolApl1, so we need to cast it to JFrame

        xs = new JSlider(SwingConstants.HORIZONTAL,0,50,10);
        xs.setMajorTickSpacing(10);
        xs.setPaintTicks(true);
            System.out.println(d3D.phase2D.atomCount);
            cf = new ChemFrame(d3D.phase2D.atomCount);
            String aname ="H";
            try{int i=0;
                for(Atom a=d3D.phase2D.firstAtom(); a!=null; a=a.nextAtom()) {
                    cf.addVert(aname, 
                        10.f*(float)a.coordinate().position(0),
                        10.f*(float)a.coordinate().position(1),
                        10.f*(float)a.coordinate().position(2));
                    i++;
 //               cf.addVert(aname, 0,0,0);
 //               cf.addVert(aname, 1,1,1);
                }
            }
            catch (Exception e)
            {
            }
                  
        myJmolPanel.display.md =cf;

        xs.addChangeListener( new ChangeListener(){
                public void stateChanged(ChangeEvent e)
                { 
                    float f= xs.getValue();
                    cf.vert[0]=f;
                    myJmolPanel.display.repaint();   
                }
            }
        );
         
        Container c=getContentPane();
        c.add(xs,BorderLayout.SOUTH);
        c.add(myJmolPanel,BorderLayout.CENTER);
        setSize(600,700);
        show();
    }  //end of constructor
  }    //end of JmolDisplay
}
