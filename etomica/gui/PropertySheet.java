package etomica.gui;

import java.beans.*;
import java.lang.reflect.*;
import java.awt.*;
import java.util.EventObject;
import java.awt.event.MouseEvent;
import java.util.Vector;
import javax.swing.DefaultCellEditor;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.JTable;
import javax.swing.tree.*;
import javax.swing.table.*;
import etomica.gui.treetable.*;
import etomica.ConstrainedPropertyEditor;
import etomica.DimensionedDoubleEditor;
import etomica.Meter;

public class PropertySheet extends javax.swing.JInternalFrame {

    static final javax.swing.border.EmptyBorder EMPTY_BORDER = new javax.swing.border.EmptyBorder(2,4,2,2);
    static final javax.swing.Icon LEAF_ICON = new javax.swing.ImageIcon("clearpixel.gif");
    static final javax.swing.Icon OPEN_ICON = new javax.swing.ImageIcon("TreeExpanded.gif");
    static final javax.swing.Icon CLOSED_ICON = new javax.swing.ImageIcon("TreeCollapsed.gif");
    
    private PropertySheetPanel panel;
    private boolean started;
    private Panel debugPanel;
    /**
     * Period (in milliseconds) between successive updates of displayed properties that might
     * change while the simulation is running (e.g., meter averages).
     * Default is 1000 msec (1 second).
     */
    public static int updateSleepPeriod = 1000;
        
    public PropertySheet(Wrapper target, int x, int y) {
	    super("Properties - <initializing...>",
	            true, //resizable
	            true  //closable
	          );
	    getContentPane().setLayout(new GridLayout());
	    
	    setBackground(Color.lightGray);	
        setBounds(x,y, 600, 300);
	    setMaximizable(true);
	    setClosable(true);
	    setIconifiable(true);
	    setVisible(true);
        
	    panel = new PropertySheetPanel(this);
	    getContentPane().add(panel);
        
	    panel.setTarget(target);
	    setTitle("Properties - " + target.getBeanLabel());

        //turn off updateThread when closing window.
        addInternalFrameListener(new javax.swing.event.InternalFrameAdapter() {
            public void internalFrameClosing(javax.swing.event.InternalFrameEvent evt) {panel.updateThread.clear();}
        });
	    started = true;
    }

    public void setTarget(Wrapper targ) {
//	    Object bean = targ.getBean();
	    String displayName = targ.getBeanLabel();
	    panel.setTarget(targ);
	    setTitle("Properties - " + displayName);
    }

    public void setCustomizer(Customizer c) {
	    panel.setCustomizer(c);
    }

    void wasModified(PropertyChangeEvent evt) {
	    panel.wasModified(evt);
    }
}//end of PropertySheet

class PropertySheetPanel extends javax.swing.JPanel {
    JScrollPane sp,sp2;
    JTree labelTree,componentTree;
    /**
     * Instance of (inner) thread class that is used to perform updates of displayed properties that might
     * change while the simulation is running (e.g., meter averages).
     */
    UpdateThread updateThread = new UpdateThread(this);

    public static class ComboCellEditor extends DefaultCellEditor {
        private JTree tree;
        
        public ComboCellEditor(JTree tree, JComboBox c){
            super(c);
            this.tree = tree;
        }      
        
        public boolean isCellEditable(EventObject e){
            boolean rv = false;
            
            if (e instanceof MouseEvent) {
                MouseEvent me = (MouseEvent)e;
                
                if (me.getClickCount() == 1) {
                    TreePath path = tree.getPathForLocation(me.getX(), me.getY());
                    DefaultMutableTreeNode node = (DefaultMutableTreeNode)path.getLastPathComponent();
                    rv = (node.getUserObject() instanceof PropertySelector);
                }
            }
            return rv;
        }
    }//end of ComboCellEditor
    
    private static class CellEditor extends javax.swing.DefaultCellEditor implements javax.swing.table.TableCellEditor{
        private JLabel labelEditor;
        private JPanel canvasEditor;
        private JTextField textEditor;
        private JComboBox comboEditor;
        private EmptyPanel emptyPanel;
        private Component c;
        
        CellEditor () {
            super(new JTextField());  //dummy argument since there is no default constructor in superclass
            setClickCountToStart(1);
        }
        
        public Component getTableCellEditorComponent(JTable table,
            Object value, boolean selected, int row, int column) {
               return (Component)value;
        }
    }//end of CellEditor
        
    private static class CellRenderer extends DefaultTableCellRenderer {
        private JLabel labelRenderer;
        private JPanel canvasRenderer;
        private JTextField textRenderer;
        private JComboBox comboRenderer;
        private EmptyPanel emptyPanel;
        private Component c;
        
        CellRenderer () {
        }

        public Component getTableCellRendererComponent(JTable table,
            Object value, boolean selected, boolean hasFocus,
            int row, int column) {
       //     super.getTableCellRendererComponent(table, value, selected, hasFocus, row, column);
            
                c = (Component)value;
                if (c instanceof JComboBox) {
                    comboRenderer = (JComboBox)c;
                    try {
                    return new StaticTextField(comboRenderer.getSelectedItem().toString());
                    } catch(NullPointerException ex) {
                        System.out.println("Null pointer exception caught");
                    }
                    return new StaticTextField("Error here!");
                }
                else {
                    return c;
                }
        }
    }//end of CellRenderer

    PropertySheetPanel(PropertySheet frame) {  //constructor
	    this.frame = frame;
	    setLayout(null);
	    setSize(600,300);

	    // Create an event adaptor.
	    adaptor = new EditedAdaptor(frame);
    }
    
    synchronized void setTarget(Wrapper targ) {

        updateThread.clear();
        updateThread = new UpdateThread(this);
	    removeAll();
	    depth = 0;
        groups = new java.util.Vector();
	    // We make the panel invisivle during the reconfiguration
	    // to try to reduce screen flicker.

	    // As a workaround for #4056424, we avoid making the panel
	    // invisible first time though, during startup.
	    if (target != null) {
	        setVisible(false);
  	    }

	    targetWrapper = targ;
	    target = targ.getBean();
	    
        //set up the treetable
        PropertyNode rootNode = makeNode(target, null, null, null);
        if(rootNode == null) return;
        PropertyModel model = new PropertyModel(rootNode);
        JTreeTable treeTable = new JTreeTable(model);
        treeTable.getTree().setRootVisible(false);
        treeTable.getTree().putClientProperty("JTree.lineStyle","Angled");
//        ((DefaultTreeCellRenderer)treeTable.getTree().getCellRenderer()).setLeafIcon(PropertySheet.LEAF_ICON);
//        ((DefaultTreeCellRenderer)treeTable.getTree().getCellRenderer()).setClosedIcon(PropertySheet.CLOSED_ICON);
//        ((DefaultTreeCellRenderer)treeTable.getTree().getCellRenderer()).setOpenIcon(PropertySheet.OPEN_ICON);
        treeTable.setDefaultRenderer(Object.class, new CellRenderer());
        treeTable.setDefaultEditor(Object.class, new CellEditor());
  //      treeTable.setShowGrid(true);
  //      treeTable.setRowHeight(30);
  //      treeTable.setGridColor(java.awt.Color.darkGray);
    
        // Put tree in a scrollable pane
        sp = new JScrollPane(treeTable);
        BorderLayout bl = new BorderLayout();
        setLayout(bl);
        add(sp, BorderLayout.CENTER);
//        add(sp);
	    frame.getContentPane().add(this);
//        frame.getContentPane().add(sp);
	    doLayout(true);

	    processEvents = true;

	    Insets ins = frame.getInsets();
    	
	    int frameWidth = getSize().width + ins.left + ins.right + 20;
	    int frameHeight = getSize().height + ins.top + ins.bottom + 20;
        
	    frame.setSize(frameWidth,frameHeight);
//	    sp.setSize(frameWidth,frameHeight);
//	    sp2.setSize(frameWidth/2,frameHeight);
//	    frame.setSize(frameWidth,frameHeight);
//	    setLocation(ins.left, ins.top);
//	    frame.getContentPane().add(this);
        setVisible(true);
    }
    
    private PropertyNode makeNode(Object object, JLabel objectLabel, Component objectView, Component objectUnitView) {

        PropertyNode root = new PropertyNode(objectLabel, objectView, objectUnitView);
        
        //See if object can have child objects in tree (cannot if it is primitive)
        if(object == null || 
            object instanceof Number || 
            object instanceof Boolean ||
            object instanceof Character) return root;
        
        //Examine all properties of object and attach child nodes to permit properties to be edited
        
        PropertyDescriptor[] properties = null;
        
        //Introspection to get array of all properties
        BeanInfo bi = null;
        try {
	        bi = Introspector.getBeanInfo(object.getClass());
	        properties = bi.getPropertyDescriptors();
	    } 
	    catch (IntrospectionException ex) {
	        error("PropertySheet: Couldn't introspect", ex);
	        return null;
	    }

        //Collect object and its properties in a PropertyGroup and add to list of all groups
        PropertyGroup group = new PropertyGroup(object, properties);
        groups.add(group);
        
        //used to store editors of values of a meter
        DimensionedDoubleEditor averageEditor = null;
        DimensionedDoubleEditor errorEditor = null;
        DimensionedDoubleEditor currentEditor = null;
        PropertyText averageView = null;
        PropertyText errorView = null;
        PropertyText currentView = null;
        
        //Loop through all properties and determine current value and find appropriate editor
	    for (int i = 0; i < properties.length; i++) {
	        // Don't display hidden or expert properties.
	        if (properties[i].isHidden() || properties[i].isExpert()) {
		        continue;
	        }

	        Object value = null;
	        Component view = null;
	        Component unitView = null;
	        JLabel label = null;

	        String name = properties[i].getDisplayName();  //Localized display name 
	        Class type = properties[i].getPropertyType();  //Type (class) of this property
	        Method getter = properties[i].getReadMethod(); //method used to read value of property in this object
	        Method setter = properties[i].getWriteMethod();//method used to set value of property
  //          System.out.println(name);
	        // Only display read/write properties.
	        if (getter == null || (setter == null && !(object instanceof Meter)) || type == Class.class) {
		        continue;
	        }
	        try {
	            //read the current value of the property
		        Object args[] = { };
		        try {value = getter.invoke(object, args);}
		        catch(NullPointerException ex) {value = null;}
	            group.values[i] = value;

                //find and instantiate the editor used to modify value of the property
	            PropertyEditor editor = null;
	            if(properties[i].isConstrained()) {
	                editor = new ConstrainedPropertyEditor();
	            }
	            else {
        	        
        	        //property is a dimensioned number
        	        if(type == Double.TYPE) {
        	            //try to get dimension from get(property)Dimension() method
        	            etomica.units.Dimension dimension = etomica.units.Dimension.introspect(object,name,bi);
        	            //try to get dimension from getDimension() method
                        if(dimension == null) dimension = etomica.units.Dimension.introspect(object,"",bi);
        	            if(dimension != null) {
        	                editor = new DimensionedDoubleEditor(dimension);
        	            }
        	        }
        	        //property is not a dimensioned number; see if its editor was set explicitly
        	        if(editor == null) { 
	                    Class pec = properties[i].getPropertyEditorClass();
		                if (pec != null) {
		                    try {
			                    editor = (PropertyEditor)pec.newInstance();
		                    } 
		                    catch (Exception ex) {}
		                }
		            }
		            //property is not a dimensioned number and was not set explicitly
		            //have editor manager look for an appropriate editor
		            if (editor == null) {
		                editor = PropertyEditorManager.findEditor(type);
		            }
		        }//done with trying to get an editor for the property
		        
	            group.editors[i] = editor;

	            // If we can't edit this component, skip it.
	            if (editor == null) {
		            // If it's a user-defined property we give a warning.
		            String getterClass = properties[i].getReadMethod().getDeclaringClass().getName();
		            if (getterClass.indexOf("java.") != 0) {
		                System.err.println("Warning: Can't find public property editor for property \""
				        + name + "\".  Skipping.");
		            }
		            continue;
	            }

                //set the editor to the current value of the property
	            editor.setValue(value);
	            
	            //add listener that causes the wasModified method to be 
	            //invoked when editor fires property-change event
	            editor.addPropertyChangeListener(adaptor);

		    // Now figure out how to display it...
		        if (editor.isPaintable() && editor.supportsCustomEditor())
        		    view = new PropertyCanvas(frame, editor);
		        else if (editor.getTags() != null)
		            view = new PropertySelector(editor);
                else if (editor.getAsText() != null) {
		            view = new PropertyText(editor);
		        }
		        else if (editor instanceof ConstrainedPropertyEditor) {
		            view = new EmptyPanel();
		        }
		        else {
		            System.err.println("Warning: Property \"" + name 
				        + "\" has non-displayable editor.  Skipping.");
		            continue;
		        }
		        if(editor instanceof DimensionedDoubleEditor) {
		            unitView = ((DimensionedDoubleEditor)editor).unitSelector();
		            if(object instanceof etomica.Meter) {
		                ((PropertyText)view).setEditable(false);
		                if(name.equals("average")) {
		                    averageEditor = (DimensionedDoubleEditor)editor;
		                    averageView = (PropertyText)view;
		                }
		                else if(name.equals("error")) {
		                    errorEditor = (DimensionedDoubleEditor)editor;
		                    errorView = (PropertyText)view;
		                }
		                else if(name.equals("mostRecent")) {
		                    currentEditor = (DimensionedDoubleEditor)editor;
		                    currentView = (PropertyText)view;
		                }
		            }
		        }
		        else {
		            unitView = new EmptyPanel();
		        }

	        } //end of try
	        catch (InvocationTargetException ex) {
		        System.err.println("Skipping property " + name + " ; exception on target: " + ex.getTargetException());
		        ex.getTargetException().printStackTrace();
		        continue;
	        } 
	        catch (Exception ex) {
		        System.err.println("Skipping property " + name + " ; exception: " + ex);
		        ex.printStackTrace();
		        continue;
	        }

	        group.labels[i] = new MyLabel(name, Label.LEFT);
	        group.views[i] = view;
	        group.unitViews[i] = unitView;
	        
	        if(group.labels[i] != null && view != null && depth < MAX_DEPTH) {
	            depth++;
	            PropertyNode child = makeNode(group.values[i],group.labels[i],view,unitView);
	            depth--;
	            if(child != null) root.add(child);
          //      root.add(new PropertyNode(labels[i],view));
            }//end if
	        
	    }//end of loop over properties
	    
	    if(object instanceof Meter) {
	        updateThread.add((Meter)object,currentEditor,averageEditor,errorEditor,
	                                       currentView,  averageView,  errorView  );
	    }
	    if(root.isLeaf() && objectView instanceof EmptyPanel) {root = null;}
	    return root;
    }//end of makeNode method

    void stretch() {
	// This gets called when a user explicitly resizes the frame.

	    Component child = null;
	    try {
	        child = (Component)frame.getContentPane().getComponent(0);
	    } catch (Exception ex) {
	        // frame has no active children;
	        return;
	    }
	    Dimension childSize = child.getSize();
	    Dimension frameSize = frame.getSize();
	    Insets ins = frame.getInsets();
	    int vpad = ins.top + ins.bottom;
	    int hpad = ins.left + ins.right;

	    // If the frame size hasn't changed, do nothing.
	    if (frameSize.width == (childSize.width + hpad) &&
		    frameSize.height == (childSize.height + vpad)) {
	        return;
	    }

	    // We treat the new frame sizes as a future maximum for our own
	    // voluntary size changes.
	    maxHeight = frameSize.height;
	    maxWidth = frameSize.width;

	    // If we've gotten smaller, force new layout.
	    if (frameSize.width < (childSize.width + hpad) ||
			    frameSize.height < (childSize.height + vpad)) {
	        // frame has shrunk in at least one dimension.
	        setTarget(targetWrapper);
	    } else {
	    // Simply resize the contents.  Note that this won't make
	    // any ScrollPane go away, that will happen on the next
	    // focus change.
 	    child.setSize(frameSize.width - hpad, frameSize.height - vpad);
    	    
	    }
    }//end of stretch

    private void doLayout(boolean doSetSize) {
//	    if (views == null || labels == null) {
//	        return;
//	    }
        if(groups.isEmpty())  return;

	    // First figure out the size of the columns.
	    int labelWidth = 92;
	    int viewWidth = 120;

	    for (java.util.Iterator iter=groups.iterator(); iter.hasNext(); ) {
	        PropertyGroup g = (PropertyGroup)iter.next();
	        for (int i = 0; i < g.labels.length; i++) {
	            if (g.labels[i] == null || g.views[i] == null) {
		        continue;
	            }
	            int w = g.labels[i].getPreferredSize().width;
	            if (w > labelWidth) {
		        labelWidth = w;
	            }
	            w = g.views[i].getPreferredSize().width;
	            if (w > viewWidth) {
		        viewWidth = w;
	            }
	        }
	    }
	    int width = 3*hPad + labelWidth + viewWidth;

	    // Now position all the components.
	    PropertyGroup g = (PropertyGroup)groups.elementAt(0);
	    int y = 10;
	    for (int i = 0; i < g.labels.length; i++) {
	        if (g.labels[i] == null || g.views[i] == null) {
		    continue;
	        }
	        g.labels[i].setBounds(hPad, y, labelWidth, 20);

//	        Dimension viewSize = views[i].getPreferredSize();
	        int h = 20;//viewSize.height;
	        //if (h < 20) {
		    //    h = 20;
	        //}
	        g.views[i].setBounds(labelWidth + 2*hPad, y, viewWidth, h);
	        y += (h + vPad);
	    }

	    y += vPad;

	    if (doSetSize) {
	        setSize(width, y);
	    }
    }

/*    public void doLayout() {
	    doLayout(false);
    }*/

    synchronized void setCustomizer(Customizer c) {
	    if (c != null) {
	        c.addPropertyChangeListener(new EditedAdaptor(frame));
	    }
    }

    /**
     * When property event is fired, this method transmits the change to the actual
     * instance of the object being edited.  Also checks for change in value of any
     * other properties and repaints the sheet accordingly.
     */
    synchronized void wasModified(PropertyChangeEvent evt) {
        if (!processEvents) {
	        return;
	    }

	    if (evt.getSource() instanceof PropertyEditor) {
	        PropertyEditor editor = (PropertyEditor) evt.getSource();
	        for (java.util.Iterator iter=groups.iterator(); iter.hasNext(); ) {
	            PropertyGroup g = (PropertyGroup)iter.next();
	            
	            for (int i = 0 ; i < g.editors.length; i++) {
	                if (g.editors[i] == editor) {
		                PropertyDescriptor property = g.properties[i];
		                Object value = editor.getValue();
		                g.values[i] = value;
		                Method setter = property.getWriteMethod();
		                if(setter == null) continue;
		                try {
		                    Object args[] = { value };
		                    args[0] = value;
		                    setter.invoke(g.target, args);
            		        
		                    // We add the changed property to the targets wrapper
		                    // so that we know precisely what bean properties have
		                    // changed for the target bean and we're able to
		                    // generate initialization statements for only those
		                    // modified properties at code generation time. 
                            targetWrapper.getChangedProperties().addElement(g.properties[i]);

		                } 
		                catch (InvocationTargetException ex) {
		                    if (ex.getTargetException() instanceof PropertyVetoException) {
			                    //warning("Vetoed; reason is: " 
			                    //        + ex.getTargetException().getMessage());
			                    // temp dealock fix...I need to remove the deadlock.
			                    System.err.println("WARNING: Vetoed; reason is: " 
					                    + ex.getTargetException().getMessage());
		                    }
		                    else error("InvocationTargetException while updating " + property.getName(), ex.getTargetException());
		                } 
		                catch (Exception ex) {
		                    error("Unexpected exception while updating " 
		                            + property.getName(), ex);
	                    }
	                    if(evt.getSource() instanceof DimensionedDoubleEditor) {
	            //            g.views[i].repaint();
	                    }
		                break;  //should rewrite this to break out of Iterator loop
		            }//end of for(int loop
		        }//end of for(Iterator loop
	        }
	    }

	    // Now re-read all the properties and update the editors
	    // for any other properties that have changed.
	    for (java.util.Iterator iter=groups.iterator(); iter.hasNext(); ) {
	        PropertyGroup g = (PropertyGroup)iter.next();
	        for (int i = 0; i < g.properties.length; i++) {
	            Object o;
	            Method setter = null;
	            try {
	                Method getter = g.properties[i].getReadMethod();
	                setter = g.properties[i].getWriteMethod();
	                Object args[] = { };
	                o = getter.invoke(g.target, args);
	            } 
	            catch (Exception ex) { o = null; }
	            if (o == g.values[i] || (o != null && o.equals(g.values[i])) ||  setter == null) {
	                // The property is equal to its old value.
    		        continue;
	            }
	            g.values[i] = o;
	            // Make sure we have an editor for this property...
	            if (g.editors[i] == null)
		            continue;
    	        
	            // The property has changed!  Update the editor.
	            g.editors[i].setValue(o);
	            if (g.views[i] != null)
		            g.views[i].repaint();
	        }//end of for(int i) loop

	        // Make sure the target bean gets repainted.
	        if (Beans.isInstanceOf(g.target, Component.class)) {
	            ((Component)(Beans.getInstanceOf(g.target, Component.class))).repaint();
	        }
	        
	        //this ensures display is updated if editor changes units
	        if(evt.getSource() instanceof DimensionedDoubleEditor) repaint();
	    }//end of for(java.util.Iterator) loop
    }//end of wasModified method

    private void warning(String s) {
    //	new ErrorDialog(frame, "Warning: " + s);
    }

    //----------------------------------------------------------------------
    // Log an error.
    private void error(String message, Throwable th) {
//	    String mess = message + ":\n" + th;
	    System.err.println(message);
	    th.printStackTrace();
	    // Popup an ErrorDialog with the given error message.
    //	new ErrorDialog(frame, mess);

    }

    private class MyLabel extends JLabel implements java.awt.event.ActionListener {
        String label;
        
        MyLabel(String l, int horizPos){
            super(l, horizPos);
            label = l;
        }
        
        public String toString(){ return label; }
        
        //causes update of sheet when unit changes for some dimensioned property
        public void actionPerformed(java.awt.event.ActionEvent evt) {
//            PropertySheetPanel.this.repaint();
        }
    }
    
    private static class PropertyGroup {
        Object target;
        PropertyDescriptor[] properties;
        PropertyEditor[] editors;
        Object[] values;
        Component[] views;
        JLabel[] labels;
        Component[] unitViews;
        
        PropertyGroup(Object t, PropertyDescriptor[] p) {
            target = t;
            properties = p;
	        editors = new PropertyEditor[properties.length];
	        values = new Object[properties.length];
	        views = new Component[properties.length];
	        labels = new JLabel[properties.length];
	        unitViews = new Component[properties.length];
	    }
    }
    
    private static class EmptyPanel extends JTextField {
        EmptyPanel() {
            super("");
            setBorder(PropertySheet.EMPTY_BORDER);
        }
    }
    
    /**
     * Used to display in the property sheet the current value of a combobox editor.
     */
    private static class StaticTextField extends JTextField {
        StaticTextField(String s) {
            super(s);
            setEditable(false);
            setBorder(PropertySheet.EMPTY_BORDER);
            setOpaque(false);
        }
    }
    
    /**
     * Class used to create a thread that updates display of properties that change
     * as the simulation progresses, and not because of user actions.  
     * These usually are values reported by meters and which are displayed on the property sheet.
     */
    static class UpdateThread extends Thread {
        private UpdateGroup first = null;
        private boolean running = false;
        private PropertySheetPanel sheet;
        UpdateThread(PropertySheetPanel s) {super(); sheet = s;}
        synchronized void add(Meter m, DimensionedDoubleEditor c, DimensionedDoubleEditor a, DimensionedDoubleEditor e,
                              PropertyText cV, PropertyText aV, PropertyText eV) {
            first = new UpdateGroup(m,c,a,e,cV,aV,eV,first);
            if(!running) {
                running = true;
                this.start();
            }
        }
        synchronized void clear() {first = null;}
        
        public void run() {
            while(first != null) {
                //linked-list loop through all meter/editor groups
                for(UpdateGroup group=first; group!=null; group=group.next) {
                    group.current.setValue(group.meter.mostRecent());
                    group.average.setValue(group.meter.average());
                    group.error.setValue(group.meter.error());
                    //update textfields values
                    group.currentView.propertyChange(null);
                    group.averageView.propertyChange(null);
                    group.errorView.propertyChange(null);
                }
                sheet.repaint();
                try { Thread.sleep(PropertySheet.updateSleepPeriod); }
                catch (InterruptedException e) { }
            }//end of while
            running = false;
        }//end of run
        
        //data structure that collects together a meter and all the editors for its values
        //also holds a handle to another UpdateGroup instance to facilite a linked list
        private static class UpdateGroup {
            Meter meter;
            UpdateGroup next;
            DimensionedDoubleEditor current, average, error;
            PropertyText currentView, averageView, errorView;
            UpdateGroup(Meter m, //constructor
                        DimensionedDoubleEditor c, 
                        DimensionedDoubleEditor a, 
                        DimensionedDoubleEditor e,
                        PropertyText cV,
                        PropertyText aV,
                        PropertyText eV,
                        UpdateGroup n) {
                meter = m;
                current = c;
                average = a;
                error = e;
                currentView = cV;
                averageView = aV;
                errorView = eV;
                next = n;
            }//end of constructor
        }//end of UpdateGroup class
        
    }//end of UpdateThread class
    
    //----------------------------------------------------------------------
    private PropertySheet frame;
    private EditedAdaptor adaptor;
    
    // We need to cache the targets' wrapper so we can annoate it with
    // information about what target properties have changed during design
    // time.
    private Wrapper targetWrapper;   
    private Object target;
    private java.util.Vector groups = new java.util.Vector();
    private boolean processEvents;
    private static int hPad = 4;
    private static int vPad = 4;
    private int maxHeight = 500;
    private int maxWidth = 500;
    
    //maximum tree depth for property nesting
    static final int MAX_DEPTH = 3;
    private int depth = 0;
    
}
