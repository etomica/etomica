package simulate;
import java.awt.*;
import java.beans.Beans;
import java.awt.event.*;

    public abstract class Display extends Canvas implements simulate.IntegrationIntervalListener, MouseListener, MouseMotionListener, KeyListener {

	Simulation parentSimulation;
	Phase phase;
    int pixels = 300;
    Image offScreen;
    Graphics osg;
    int updateInterval;
    int iieCount;
    Component displayTool = null;  //displayTool is some component inside a Display object that is doing all the painting (not often used)
    private Display nextDisplay;
    private Display previousDisplay;
    
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
	
	int offset = 10;  //width of region at perimeter at which move or resize is indicated
	boolean moving = false;
	boolean resizing = false;
	boolean mouseInside = false;
	int cursorLocation = DEFAULT_AREA;
	Point anchor;
	Point anchorAbs= new Point();
	Point initialLocation;
	Dimension initialSize;
	
    boolean resizable = false;  //flag to indicate if display can be resized
    boolean movable = false;     //flag to indicate if display can be moved
    
    public Display () {
        setSize(pixels, pixels);
        setBackground(Color.white);
	    setUpdateInterval(1);
        addMouseListener(this);
        addMouseMotionListener(this);
        addKeyListener(this);
    }
    
    public void setMovable(boolean b) {movable = b;}
    public boolean isMovable() {return movable;}
    public void setResizable(boolean b) {resizable = b;}
    public boolean isResizable() {return resizable;}

    public final Display getNextDisplay() {return nextDisplay;}
    public final Display getPreviousDisplay() {return previousDisplay;}
   /**
    * Sets the display following this one in the linked list of displays.
    *
    * @param d the display to be designated as this display's nextDisplay
    */
    public final void setNextDisplay(Display d) {
      this.nextDisplay = d;
      d.previousDisplay = this;
    }
    
    public void setPhase(Phase p) {phase = p;}
    
    public void createOffScreen () {
        if (offScreen == null) { 
            createOffScreen(pixels);
        }
    }
    public void createOffScreen (int p) {
        createOffScreen(p,p);
    }
    public void createOffScreen(int w, int h) {
        offScreen = createImage(w,h);
        osg = offScreen.getGraphics();
    }

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

    public void integrationIntervalAction(IntegrationIntervalEvent evt) {
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        doUpdate();
            repaint();
	    }
    }

    public abstract void doUpdate();
    
    public abstract void doPaint(Graphics g);
    
	public void setParentSimulation(Simulation s) {parentSimulation = s;}

    public final int getUpdateInterval() {return updateInterval;}
    public final void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
    }
    
    // MouseListener methods
    
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

    //keyListener methods
    
    public void keyPressed(KeyEvent e) {System.out.println("keypressed:");}
    public void keyReleased(KeyEvent e) {System.out.println("keyreleased: ");}
    public void keyTyped(KeyEvent e) {
        char c = e.getKeyChar();
        System.out.println("keytyped: "+c);
        if(Character.isDigit(c)) {digitTyped(Character.getNumericValue(c));}
        else if(Character.isLetter(c)) {letterTyped(c);}
    }
    
    protected void digitTyped(int i) {}   //override if want to process number keypress
    protected void letterTyped(char c) {} //override if want to process letter keypress
}
