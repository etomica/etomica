
package etomica.gui;

/**
 * This wrapper class keeps track of various BeanBox related
 * state for each bean in the composition window.
 *
 * Among other things, it draws the black-and-white hashed
 * border around the currenty active component.
 */

import java.awt.*;
import java.awt.Component;
import java.awt.event.*;
import java.io.*;
import java.beans.*;
import java.util.*;
import java.lang.reflect.Method;
import java.lang.reflect.InvocationTargetException;
import java.applet.Applet;
import javax.swing.JComponent;
import javax.swing.JPanel;
import sunw.beanbox.PropertyHookup;

public class Wrapper extends JPanel implements Serializable {
    
    static final long serialVersionUID = 1144602051002987355L;

    public Wrapper(Object bean, String beanLabel, String beanName) {
        this.bean = bean;
	    if (beanName == null) {
	        beanName = bean.getClass().getName();
	    }
	    if (beanLabel == null) {
	        beanLabel = beanName;
	    }
	    this.beanName = beanName;
	    this.beanLabel = beanLabel;
	    this.isFromPrototype = false;
    	
	    setLayout(null);
	    if (Beans.isInstanceOf(bean, JComponent.class)) {
            child = (JComponent)Beans.getInstanceOf(bean, JComponent.class);
	    } 
	    else {
	        invisibleWrappers.addElement(this);
	        child = new javax.swing.JLabel(beanLabel);
	    }
	    add(child);
	    child.setLocation(borderWidth, borderWidth);
	    initialize();
	    attachListeners();
    }

    void initialize() { 
	    if (! Beans.isInstanceOf(bean, JComponent.class)) {
	        invisibleWrappers.addElement(this);
	    }
	    esdMap = new Hashtable();
	    try {
	        BeanInfo bi = Introspector.getBeanInfo(bean.getClass());
	        EventSetDescriptor[] esds = bi.getEventSetDescriptors();
	        for (int i = 0; i < esds.length; i++) {
	            esdMap.put(esds[i].getName(), esds[i]);
	        }
	    } 
	    catch (IntrospectionException ex) {
	        System.err.println("Wrapper couldn't introspect on bean: " + bean);
        }
	    if (eventTargets == null) {
	        eventTargets = new Vector();
	        propertyTargets = new Vector();
	    } 
	    else {
	        propertyTargets = getWPEIfromWET(eventTargets);
	    }
    }

    /**
     * This method massages information already present in WrapperEventTarget
     * into Wrapper(Property)EventInfo, which is easier to use.
     */

    private Vector getWPEIfromWET(Vector wets) {
	    Vector back = new Vector();
	    // we only care about properties, not method hookups
	    for (Enumeration e = wets.elements(); e.hasMoreElements(); ) {
            WrapperEventTarget wet = (WrapperEventTarget) e.nextElement();
	        
	        if (wet.targetListener instanceof PropertyHookup) {
		        PropertyHookup h = (PropertyHookup) wet.targetListener;
                Hashtable table = h.getTargetsByProperty();
		        for (Enumeration keys = table.keys(); keys.hasMoreElements(); ) {
		            String propertyName = (String) keys.nextElement();
		            Vector targets = (Vector) table.get(propertyName);

		            for (Enumeration ee = targets.elements(); ee.hasMoreElements(); ) {
			            Object t = ee.nextElement();
			            back.addElement( new WrapperPropertyEventInfo(h.getTargetObject(t), propertyName, h.getSetterMethod(t)));
		            }
		        }
	        }
	    }
	    return back;
    }


    /**
     * Serialization methods
     */

    private void readObject(ObjectInputStream s) throws ClassNotFoundException, IOException {
	    s.defaultReadObject();
	    initialize();
    }

    synchronized void addPropertychangeListener(String propertyName, PropertyChangeListener listener) {}

    /**
     * How many Hookups?
     */
    int getEventHookupCount() {
	    if (propertyTargets.size() > 0) {
	        return eventTargets.size() + propertyTargets.size() - 1;
	    } 
	    else {
	        return eventTargets.size();
	    }
    }

    public String getAdderName(String eventSetName) {
	    EventSetDescriptor esd = (EventSetDescriptor) esdMap.get(eventSetName);
	    Method adder;
	    adder = esd.getAddListenerMethod();
	    return adder.getName();
    }

    public String getRemoverName(String eventSetName) {
	    EventSetDescriptor esd = (EventSetDescriptor) esdMap.get(eventSetName);
	    Method remover;
	    remover = esd.getRemoveListenerMethod();
	    return remover.getName();
    }

    /**
     * This will replace the bottom two
     */
    public WrapperEventInfo[] getEventHookupInfo() {
    	int s = getEventHookupCount();

	    WrapperEventInfo[] back = new WrapperEventInfo[s];
	    int i = 0;
	    Enumeration e = eventTargets.elements();
        
	    while (e.hasMoreElements()) {
	        WrapperEventTarget et = (WrapperEventTarget) e.nextElement();
	        if (et.targetBean == null) {
		        // this is a property bound hookup
	        } 
	        else {
		        back[i] = new WrapperEventInfo(et.targetBean, et.targetListener.getClass().getName(), et.eventSetName);
		        i += 1;
	        }
	    }
	    e = propertyTargets.elements();
	    while (e.hasMoreElements()) {
	        back[i] = (WrapperEventInfo) e.nextElement();
	        i += 1;
	    }
	    return back;
    }


    // Add a (set of) PropertyHookup
    synchronized void addPropertyTarget(String propertyName, Object targetObject, Method setter) {
	    propertyTargets.addElement( new WrapperPropertyEventInfo(targetObject, propertyName, setter));
    }

    // Add a hookup.  All property bound hookups are represented by (at most) one hookup
    synchronized void addEventTarget(String eventSetName, Wrapper targetWrapper, Object listener) {
	    WrapperEventTarget et = new WrapperEventTarget();
	    et.eventSetName = eventSetName;
	    if (targetWrapper != null) {
	        et.targetBean = targetWrapper.getBean();
	    }
	    et.targetListener = listener;
	    eventTargets.addElement(et);
	    EventSetDescriptor esd = (EventSetDescriptor) esdMap.get(eventSetName);
	    if (esd == null) {
	        System.err.println("Internal error: Wrapper.addEventTarget missing event set");
	        System.err.println("        eventSetName = " + eventSetName);
	        System.err.println("        bean = " + bean);
	        return;
	    }
	    Method adder = esd.getAddListenerMethod();
	    Method remover = esd.getRemoveListenerMethod();
	    if (adder == null || remover == null) {
	        System.err.println("Internal error: Wrapper.addEventTarget missing add/remote listener");
	        System.err.println("        eventSetName = " + eventSetName);
	        System.err.println("        bean = " + bean);
	        return;
	    }
	    try {
	        Object args[] = { listener };
	        adder.invoke(bean, args);
	    } 
	    catch (InvocationTargetException ex) {
	        System.err.println("Wrapper: adding event listener for " + eventSetName + " failed:");
	        System.err.println("    " + ex.getTargetException());
	    } 
	    catch (Exception ex) {
	        System.err.println("Wrapper: adding event listener for " + eventSetName + " failed:");
	        System.err.println("    " + ex);
	    }
    }

    /**
     * Temporarily remove any event listeners.
     */

    void removeListeners() {
	    Enumeration enum = eventTargets.elements();
	    while (enum.hasMoreElements()) {
	        WrapperEventTarget et = (WrapperEventTarget)enum.nextElement();
	        EventSetDescriptor esd = (EventSetDescriptor) esdMap.get(et.eventSetName);
	        Method remover = esd.getRemoveListenerMethod();
	        try {
	            Object args[] = { et.targetListener };
	            remover.invoke(bean, args);
	        } 
	        catch (InvocationTargetException ex) {
	            System.err.println("Wrapper: removing event listener for "
					    + et.eventSetName + " failed:");
	            System.err.println("    " + ex.getTargetException());
		        ex.getTargetException().printStackTrace();
	        } 
	        catch (Exception ex) {
	            System.err.println("Wrapper: removing event listener for "
					    + et.eventSetName + " failed:");
	            System.err.println("    " + ex);
		        ex.printStackTrace();
	        }
	    }
	    // Remove mouse listeners.
	    listenForMice(false);
    }

    /**
      * (Re)-attach any event listeners.
      */
    void attachListeners() {   
	    Enumeration enum = eventTargets.elements();
	    while (enum.hasMoreElements()) {
	        WrapperEventTarget et = (WrapperEventTarget)enum.nextElement();
	        EventSetDescriptor esd = (EventSetDescriptor) esdMap.get(et.eventSetName);
	        Method adder = esd.getAddListenerMethod();
	        try {
	            Object args[] = { et.targetListener };
	            adder.invoke(bean, args);
	        } 
	        catch (InvocationTargetException ex) {
	            System.err.println("Wrapper: adding event listener for "
					    + et.eventSetName + " failed:");
	            System.err.println("    bean = " + bean);
	            System.err.println("    " + ex.getTargetException());
		        ex.getTargetException().printStackTrace();
	        } 
	        catch (Exception ex) {
	            System.err.println("Wrapper: adding event listener for "
					    + et.eventSetName + " failed:");
	            System.err.println("    bean = " + bean);
	            System.err.println("    " + ex);
		        ex.printStackTrace();
	        }
	    }
	    // Reattach mouse listeners.
	    listenForMice(true);
    }


    void writeNakedBean(ObjectOutputStream oos) throws IOException {
	    // First, detach all event listeners.
	    removeListeners();
	    try { //Now write the bean.
	        oos.writeObject(bean);
	    } 
	    finally { // Now rettach all event listeners.
	        attachListeners();
	    }
    }

    /**
     * Cleanup is called when a Wrapper has been "cut" from the BeanBox.
     */
    void cleanup() {
	    if (bean instanceof Applet) {
	        Applet apl = (Applet)bean;
	        apl.stop();
	        apl.destroy();
	    }
	    removeListeners();	
	    // We should also remove ourself from any event sources...
    }

    /**
     * set whether or not the target bean is from a serialized origin
     */
    public void setFromPrototype(boolean b){ isFromPrototype = b; }
	 
    /**
     * get whether or not the target bean is from a serialized origin
     */
    public boolean isFromPrototype(){ return isFromPrototype; }

    /**
     * get the wrapped bean.
     */
    public Object getBean() { return bean; }

    /**
     * get the wrapped beanName
     */
    public String getBeanLabel() { return beanLabel; }

    /**
     * get the wrapped beanName
     */
    public String getBeanName() { return beanName; }

    /**
     * get the AWT component used to represent the wrapped bean.
     */
    public JComponent getChild() { return child; }

    /**
     * get the properties changed at design time.
     */
    public Vector getChangedProperties()
    {
	    if (changedProperties == null) 
	        changedProperties = new Vector();
	    return changedProperties;
    }
	 
    public void doLayout() {
	    // Has the child gotten bigger?  If so, expand to fit.
	    Dimension d = getSize();
	    Dimension cd = child.getMinimumSize();
	    if (cd.width > (d.width - (2*borderWidth)) ||
		    cd.height > (d.height - (2*borderWidth))) {
	        int width = d.width;
	        if (cd.width > (d.width - (2*borderWidth)))
		        width = cd.width + (2*borderWidth);
	        int height = d.height;
	        if (cd.height > (d.height - (2*borderWidth))) 
	            height = cd.height + (2*borderWidth);
	        setSize(width,height);
	    }
    }

/*    public void setActive(boolean isActive) {
	    active = isActive;
	    repaint();
    }
*/
    public Dimension getPreferredSize() {
	    Dimension childSize = child.getPreferredSize();
	    if (childHasStupidPreferredSize()) {
	        childSize = child.getSize();
	    }
	    int width = childSize.width;
	    int height = childSize.height;
	    // Make sure the child is at least its minimum size.
	    Dimension minSize = child.getMinimumSize();
	    if (minSize.height > height) {
	        height = minSize.height;
	    }
	    if (minSize.width > width) {
	        width = minSize.width;
	    }
	    width += (2 * borderWidth);
	    height += (2 * borderWidth);
	    return new Dimension(width, height);
    }

    private boolean childHasStupidPreferredSize() {
	    // We do a special check for demented behaviour from empty Panels
	    // for the sake of applets.  If an applet does not explicitly
	    // specify a dimension, then we allow it to be resized arbitrarily
	    // rather than limiting it to the tiny size returned by the FlowLayout
	    if (child instanceof JPanel) {
	        Container cont = (Container)child;
	        LayoutManager lay = cont.getLayout();
	        if (cont.getComponentCount() == 0 && lay instanceof FlowLayout) {
		        FlowLayout flow = (FlowLayout) lay;
		        Dimension cd = child.getPreferredSize();
		        Dimension fd = flow.preferredLayoutSize(cont);
		        if (cd.width == fd.width && cd.height == fd.height) {
		            return true;
		        }
	        }
	    }
	    return false;
    }

    public void setBounds(int x, int y, int width, int height) {

	    // If we're the top level wrapper, there is no dickering.
	    if (getParent() != null && getParent() instanceof Frame) {
	        super.setBounds(x,y,width,height);
	        child.setBounds(borderWidth, borderWidth,
		    width-(2*borderWidth), height-(2*borderWidth));
	        child.validate();
	        return;
	    }

	    // Figure out what size we want to set the child 
	    width -= (2 * borderWidth);
	    height -= (2 * borderWidth);

	    // Make sure the child is at least its minimum size.
	    Dimension minSize = child.getMinimumSize();
	    if (minSize.height > height) {
	        height = minSize.height;
	    }
	    if (minSize.width > width) {
	        width = minSize.width;
	    }

	    // Make sure the child is under its maximum size.
	    Dimension maxSize = child.getMaximumSize();
	    if (height > maxSize.height) {
	        height = maxSize.height;
	    }
	    if (width > maxSize.width) {
	        width = maxSize.width;
	    }

	    // Now we can set the child's size.
	    child.setBounds(borderWidth, borderWidth, width, height);

	    // Finally we can set our own size.
	    width += (2 * borderWidth);
	    height += (2 * borderWidth);
	    super.setBounds(x, y, width, height);

	    child.validate();
    }

    // Note that to avoid deadlocks, paint is *not* synchronized
    public void paint(Graphics g) {
	    if (active && Beans.isDesignTime()) {
	        int width = getSize().width;
	        int height = getSize().height;

	        getHashBars(this);

	        // Draw the bounding hasbox as a set of images.

	        // First draw the top and bottom bars.
	        int nudge = 2 * hashBarWidth;
	        int bottomNudge = - ((height-hashBarWidth) % nudge);
	        for (int x = 0; x < width; x += hashBarLength) {
	            g.drawImage(xHashBar, x, 0, null);
	            g.drawImage(xHashBar, x + bottomNudge, height-hashBarWidth, null);
	        }
	        // Now draw the left and right bars.
	        int rightNudge = - ((width-hashBarWidth) % nudge);
	        for (int y = 0; y < height; y += hashBarLength) {
	            g.drawImage(yHashBar, 0, y, null);
	            g.drawImage(yHashBar, width-hashBarWidth, y+rightNudge, null);
	        }


	    }
	    super.paint(g);
    }

    //----------------------------------------------------------------------
    private static Hashtable eventModelCache = new Hashtable();

    // Check if our child needs to use the old AWT event mnodel.

    private synchronized boolean useNewEventModel() {
	    // Check our cache first.
	    Boolean b = (Boolean) eventModelCache.get(child.getClass());
	    if (b != null) {
	        return b.booleanValue();
	    }

	    // We check whether the bean has any public methods that take
	    // old AWT event objects.  If so, we assume it wants the
            // old event model, otherwise we assume it wants the new.
	    boolean useNew = true;

	    try {
	        Class clz = child.getClass();
	        while (useNew && clz != null) {

		    java.lang.reflect.Method methods[] = clz.getDeclaredMethods();

		    for (int i = 0; i < methods.length; i++) {
		        java.lang.reflect.Method m = methods[i];
		        int mods = m.getModifiers();
		        if (!java.lang.reflect.Modifier.isPublic(mods)) {
			    // Skip non=public methods.
			    continue;
		        }
		        Class params[] = m.getParameterTypes();
		        if (params.length > 0 && params[0] == java.awt.Event.class) {
			    // First arg is java.awt.Event.  We assume this is an
			    // old-style AWT event handler method.
			    useNew = false;
			    break;
		        }
		    }
		    clz = clz.getSuperclass();
		    if (clz.getName().indexOf("java.") == 0) {
		    // We've reached a java.* class, so we're done.
		        break;
	  	    }
	        }
	    } catch (Exception ex) {
	        System.err.println("Wrapper.useOldEventModel caught: " + ex);
	    }

	    eventModelCache.put(child.getClass(), new Boolean(useNew));
	    return useNew;
    }


    void listenForMice(boolean enable) {
	    // If the child uses the old event mode we rely on receiving events
	    // through handleEvent.  Otherwise we use the new event model and
	    // register explicit listeners for mouse events.
	    if (!useNewEventModel()) {
	        System.err.println("WARNING: \"" 
							    + child.getClass().getName()
							    + "\" "
							    + "is a transitional bean.\n"
							    + "SOME BEAN CONTAINERS MAY NOT SUPPORT"
							    + " TRANSITIONAL BEANS!"
							    );
	        return;
	    }
    }

    private static synchronized void getHashBars(JComponent c) {
	    if (xHashBar != null) {
	        return;
	    }
	    int len = hashBarLength + 20;

	    xHashBar = c.createImage(len, hashBarWidth);
	    yHashBar = c.createImage(hashBarWidth, len);

	    Polygon poly = new Polygon();
	    Graphics g = xHashBar.getGraphics();
	    for (int i = 0; i < 4; i++) {
	        poly.addPoint(0,0);
	    }
   	    poly.ypoints[2] = hashBarWidth;
   	    poly.ypoints[3] = hashBarWidth;

	    for (int x = 0; x < (hashBarWidth+len); x += hashBarWidth) {
	        // draw alternate dark gray and light gray stripes.
	        if (((x/hashBarWidth)%2) == 0) {
	            g.setColor(Color.darkGray);
	        } else {
	            g.setColor(Color.lightGray);
	        }
	        poly.xpoints[0] = x;
	        poly.xpoints[1] = x + hashBarWidth;
	        poly.xpoints[2] = x;
	        poly.xpoints[3] = x - hashBarWidth;
	        g.fillPolygon(poly);
	    }

	    g = yHashBar.getGraphics();
   	    poly.xpoints[0] = 0;
   	    poly.xpoints[1] = hashBarWidth;
   	    poly.xpoints[2] = hashBarWidth;
   	    poly.xpoints[3] = 0;

	    for (int y = 0; y < (hashBarWidth+len); y += hashBarWidth) {
	        // draw alternate dark gray and light gray stripes.
	        if (((y/hashBarWidth)%2) == 0) {
	            g.setColor(Color.darkGray);
	        } else {
	            g.setColor(Color.lightGray);
	        }
	        poly.ypoints[0] = y;
	        poly.ypoints[1] = y - hashBarWidth;
	        poly.ypoints[2] = y;
	        poly.ypoints[3] = y + hashBarWidth;
	        g.fillPolygon(poly);
	    }
    }

    private static boolean isJDK(String version) { return System.getProperty("java.version").equals(version); }

    //----------------------------------------------------------------------

    private JComponent child;
    private Object bean;
    private String beanLabel;
    private String beanName;
    private transient boolean active;
    private transient Cursor cursor = defaultCursor;

    private final static int borderWidth = 5;
    private final static int resizeDelta = 8;    

    // Shorthands for the cursors.
    private static Cursor nwResizeCursor = Cursor.getPredefinedCursor(Cursor.NW_RESIZE_CURSOR);
    private static Cursor neResizeCursor = Cursor.getPredefinedCursor(Cursor.NE_RESIZE_CURSOR);
    private static Cursor swResizeCursor = Cursor.getPredefinedCursor(Cursor.SW_RESIZE_CURSOR);
    private static Cursor seResizeCursor = Cursor.getPredefinedCursor(Cursor.SE_RESIZE_CURSOR);
    private static Cursor moveCursor     = Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR);
    private static Cursor defaultCursor  = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);

    private static Vector invisibleWrappers = new Vector();
    private transient boolean sawMouseDown;

    // Table that maps events set names to EventDescriptors.
    private transient Hashtable esdMap;
    // List of WrapperEventTargets that we fire events at.
    private Vector eventTargets;
    // List of properties that we bound into
    transient private Vector propertyTargets;

    private static int hashBarWidth = 4;
    private static int hashBarLength = 104;
    private static Image xHashBar;
    private static Image yHashBar;

    // List of this wrappers' beans' changed properties
    private transient Vector changedProperties;

    private transient boolean isFromPrototype;

    //----------------------------------------------------------------------

}

// Class to hold state on event listener hookup for which this
// Wrapper's bean is a source.

class WrapperEventTarget implements Serializable {
    static final long serialVersionUID = 4831901854891942741L;

    String eventSetName;
    Object targetBean;
    Object targetListener;
}
