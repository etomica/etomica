package simulate;
import java.awt.*;
import java.beans.Beans;
import java.awt.event.*;
//import javax.swing.JComponent;

/**
 * Superclass for classes that display information from simulation by painting to a canvas.
 * Defines methods useful for dealing with mouse and key events targeted at the display.
 * Much of the class is involved with defining event handling methods to permit display 
 * to be moved or resized; in the future these functions will be handled instead using awt component functions.
 * 
 * @see DisplayPhase.Canvas
 */
public abstract class DisplayCanvas extends javax.swing.JPanel implements MouseListener, MouseMotionListener, KeyListener, java.io.Serializable {

    transient Image offScreen;
    transient Graphics osg;
    int pixels = 300;
        
    static final int DEFAULT_AREA = 0;
    static final int MOVE_AREA = 1;
    static final int RESIZE_AREA_NW = 2;
    static final int RESIZE_AREA_NE = 3;
    static final int RESIZE_AREA_SW = 4;
    static final int RESIZE_AREA_SE = 5;

	static final Cursor curDefault = new Cursor(Cursor.DEFAULT_CURSOR);
	static final Cursor curMove    = new Cursor(Cursor.MOVE_CURSOR);
	static final Cursor curResize1 = new Cursor(Cursor.NW_RESIZE_CURSOR);
	static final Cursor curResize2 = new Cursor(Cursor.NE_RESIZE_CURSOR);
    
    /**
     * Width of region at perimeter within which move or resize is indicated
     */
	int offset = 10;
	boolean moving = false;
	boolean resizing = false;
	boolean mouseInside = false;
	int cursorLocation = DEFAULT_AREA;
	Point anchor;
	Point anchorAbs= new Point();
	Point initialLocation;
	Dimension initialSize;
    
    /**
     * Flag to indicate if display can be resized
     */
    boolean resizable = false;
    /**
     * Flag to indicate if display can be moved
     */
    boolean movable = false;

    public DisplayCanvas() {
        setSize(pixels, pixels);
        setBackground(Color.white);
//        addMouseListener(this);
//        addMouseMotionListener(this);
//        addKeyListener(this);
        setMinimumSize(new Dimension(pixels, pixels));
        setMaximumSize(new Dimension(pixels, pixels));
        setPreferredSize(new Dimension(pixels, pixels));
    }
    public void createOffScreen () {
        if (offScreen == null) { 
            createOffScreen(getSize().width, getSize().height);
        }
    }
    public void createOffScreen (int p) {
        createOffScreen(p,p);
    }
    public void createOffScreen(int w, int h) {
        offScreen = createImage(w,h);
        if(offScreen != null) osg = offScreen.getGraphics();
    }
        
    public abstract void doPaint(Graphics g);
    
    public void update(Graphics g) {paint(g);}
        
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

    public void setMovable(boolean b) {movable = b;}
    public boolean isMovable() {return movable;}
    public void setResizable(boolean b) {resizable = b;}
    public boolean isResizable() {return resizable;}

    public void changeSize(int w, int h, MouseEvent e) {
        if(e.isShiftDown()) {
            int x = Math.max(w,h);
            setSize(x,x);
            createOffScreen(x);
        }
        else {
            setSize(w,h);
            createOffScreen(w,h);
        }
    }
        
    // MouseListener methods
    // Move and resize features will likely be eliminated in favor of functionality of a window pane
        
    public void mouseReleased(MouseEvent e)
    {
	    Point p=new Point(e.getX(), e.getY());
	    moving = false;
	    resizing = false;
	    if (((p.x-anchor.x)*(p.x-anchor.x)) <=4 &&
	        ((p.y-anchor.y)*(p.y-anchor.y)) <=4 )
	        {;}  //take action as if Press/Release were a click
      
    }
    public void mousePressed(MouseEvent e)
    {
	    anchor = new Point(e.getX(), e.getY());   //location of mouse when button pressed
	    moving = false;
	    resizing = false;
        if (cursorLocation != DEFAULT_AREA) {
	        if (cursorLocation == MOVE_AREA) {moving = true;}
	        else  {resizing = true;}
	        initialLocation = getLocation();
	        initialSize = getSize();
	        anchorAbs.x = anchor.x + initialLocation.x;
	        anchorAbs.y = anchor.y + initialLocation.y;
	    }
    }
        
    // Determines the location of the cursor with respect to the perimeter of the display
    // Changes cursor and prepares for action if located within offset shell about perimeter
    public void mouseMoved(MouseEvent e)
    {
        if(!movable && !resizable) return;   //no action if stationary and rigid
            
	    Rectangle r = this.getBounds();
	    int x = e.getX();
	    int y = e.getY();

	    if (x >= offset && x <= r.width-offset && y >= offset && y <= r.height-offset)  
	        {cursorLocation = DEFAULT_AREA;}
	    else if (x >= 0 && x <= offset && y >= 0 && y <= offset) 
	        {cursorLocation = RESIZE_AREA_NW;}
	    else if (x >= 0 && x <= offset && y >= r.height-offset && y <= r.height) 
	        {cursorLocation = RESIZE_AREA_SW;}
	    else if (x >= r.width-offset && x <= r.width && y >= 0 && y <= offset) 
	        {cursorLocation = RESIZE_AREA_NE;}
	    else if (x >= r.width-offset && x <= r.width  && y >= r.height-offset && y <= r.height) 
	        {cursorLocation = RESIZE_AREA_SE;}
	    else if (x >= 0 && x <= 0 +r.width && y >= 0 && y <= 0 +r.height) 
	        {cursorLocation = MOVE_AREA;}
	    else
	        {cursorLocation = DEFAULT_AREA;}
    	          
	    if (movable && cursorLocation==MOVE_AREA) setCursor(curMove);
	    else if (resizable && (cursorLocation==RESIZE_AREA_NW || cursorLocation==RESIZE_AREA_SE)) setCursor(curResize1);
	    else if (resizable && (cursorLocation==RESIZE_AREA_NE || cursorLocation==RESIZE_AREA_SW)) setCursor(curResize2);
	    else setCursor(curDefault);
            
    }
        
    public void mouseExited(MouseEvent e)  {
        transferFocus();
        mouseInside = false;
    }
    public void mouseEntered(MouseEvent e) {
        requestFocus();   //can accept key event only if this has focus
        mouseInside = true;
    }
        
    public void mouseDragged(MouseEvent e)
    {
        int deltaX = (getLocation().x + e.getX()) - anchorAbs.x;  //total x-distance moved since mousePressed
        int deltaY = (getLocation().y + e.getY()) - anchorAbs.y;  //total y-distance moved since mousePressed
	    if (moving) {
	        int ix = initialLocation.x + deltaX;
	        int iy = initialLocation.y + deltaY;
            setLocation(ix,iy);
        }
        else if (resizing) {
	        if (cursorLocation == RESIZE_AREA_NW) {deltaX *= -1; deltaY *= -1;}
	        else if (cursorLocation == RESIZE_AREA_SW) {deltaX *= -1;}
	        else if (cursorLocation == RESIZE_AREA_NE) {deltaY *= -1;}
	        else if (cursorLocation == RESIZE_AREA_SE) {;}
	        int w = initialSize.width + deltaX;
	        int h = initialSize.height + deltaY;
            changeSize(w,h,e);
        }
    }
    public void mouseClicked(MouseEvent e) {}

   /**
    * KeyListener method.  Performs no action.
    */
    public void keyPressed(KeyEvent e) {}
   /**
    * KeyListener method.  Performs no action.
    */
    public void keyReleased(KeyEvent e) {}
   /**
    * KeyListener method.  Determines if key pressed is a letter or digit, and calls digitTyped or letterTyped method, as appropriate.
    */
    public void keyTyped(KeyEvent e) {
        char c = e.getKeyChar();
        if(Character.isDigit(c)) {digitTyped(Character.getNumericValue(c));}
        else if(Character.isLetter(c)) {letterTyped(c);}
    }
    
    /**
     * Method invoked if a digit is typed while the cursor is over the canvas.
     * Default action does nothing, but may be overriden in subclasses to process key press
     * 
     * @param i is the digit typed, in integer form
     */
    protected void digitTyped(int i) {}   //override if want to process number keypress
    /**
     * Method invoked if a letter is typed while the cursor is over the canvas.
     * Default action does nothing, but may be overriden in subclasses to process key press
     * 
     * @param c is the letter typed, in character form
     */
    protected void letterTyped(char c) {} //override if want to process letter keypress
        
} //end of DisplayCanvas class

