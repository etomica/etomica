package etomica.gui;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.beans.*;
import java.lang.reflect.*;
import java.net.URL;
import java.util.EventObject;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.table.*;
import etomica.gui.treetable.*;
import etomica.ConstrainedPropertyEditor;
import etomica.DimensionedDoubleEditor;
import etomica.TypedConstantEditor;
import etomica.Meter;
import etomica.Simulation;
import etomica.units.Unit;
import javax.swing.event.TreeExpansionListener;
import javax.swing.event.TreeExpansionEvent;


public class PropertySheet extends JInternalFrame {
    
    public String getVersion() {return "PropertySheet:01.05.14";}

    static final javax.swing.border.EmptyBorder EMPTY_BORDER = new javax.swing.border.EmptyBorder(2,4,2,2);
    static ImageIcon LEAF_ICON;
    static ImageIcon OPEN_ICON;
    static ImageIcon CLOSED_ICON;
    static {
        try {
    	    LEAF_ICON = new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"clearpixel.gif"));
		    CLOSED_ICON = new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"TreeCollapsed.gif"));
		    OPEN_ICON = new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"TreeExpanded.gif"));
		}
		catch (java.net.MalformedURLException error) { }
        PropertyEditorManager.registerEditor(etomica.DataSource.class, etomica.DataSourceEditor.class);
        PropertyEditorManager.registerEditor(etomica.MCMove[].class, etomica.gui.McMoveEditor.class);
        PropertyEditorManager.registerEditor(etomica.units.Unit.class, etomica.UnitEditor.class);
		PropertyEditorManager.registerEditor(etomica.ModulatorAbstract.class, etomica.ModulatorEditor.class);
    }
    private PropertySheetPanel panel;
    private boolean started;
    private Panel debugPanel;
    /**
     * Period (in milliseconds) between successive updates of displayed properties that might
     * change while the simulation is running (e.g., meter averages).
     * Default is 1000 msec (1 second).
     */
    public static int updateSleepPeriod = 1000;
        
    public PropertySheet(Simulation.Element element, int x, int y) {
	    super("Properties - <initializing...>",
	            true, //resizable
	            true  //closable
	          );
	    getContentPane().setLayout(new GridLayout());
	    
	    setBackground(Color.lightGray);	
        setBounds(x,y, 400, 300);
	    setMaximizable(true);
	    setClosable(true);
	    setIconifiable(true);
	    setVisible(true);
        
	    panel = new PropertySheetPanel(this);
	    getContentPane().add(panel);
        
        setTarget(element);
        
        //turn off updateThread when closing window.
        addInternalFrameListener(new javax.swing.event.InternalFrameAdapter() {
            public void internalFrameClosing(javax.swing.event.InternalFrameEvent evt) {
                panel.updateThread.clear();
            }
        });
	    started = true;
    }

    /**
     * Sets the property sheet for editing the properties of the given simulation element.
     */
    public void setTarget(Simulation.Element element) {
        setSize(250,300);///
	    panel.setTarget(element);
        if (element != null && element.parentSimulation() != null)
	        setTitle("Properties - " + element.parentSimulation().getName());
        else if(Simulation.instance != null) setTitle("Properties - " + Simulation.instance.getName());
        else setTitle("Properties - " + "null");
    }

    void wasModified(PropertyChangeEvent evt) {
	    panel.wasModified(evt);
    }
}//end of PropertySheet


/**
 * This is the panel forming the property sheet, and which is placed inside the property-sheet frame.
 * It displays a ScrollPane which houses a JTreeTable.
 */
class PropertySheetPanel extends JPanel {
    JScrollPane sp,sp2;
    JTreeTable treeTable;
    PropertyModel model;
    /**
     * Instance of (inner) thread class that is used to perform updates of displayed properties that might
     * change while the simulation is running (e.g., meter averages).
     */
    UpdateThread updateThread = new UpdateThread(this);


    /**
     * Constructs a property sheet panel for placement into the given frame.
     */
    PropertySheetPanel(PropertySheet frame) {  //constructor
	    this.frame = frame;
	    setLayout(null);
	    setSize(250,300);

	    // Create an event adaptor.
	    adaptor = new EditedAdaptor(frame);
    }
    
    private static class CellEditor extends DefaultCellEditor implements javax.swing.table.TableCellEditor{
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
        
        CellRenderer () {}

        public Component getTableCellRendererComponent(JTable table,
            Object value, boolean selected, boolean hasFocus,
            int row, int column) {
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
                else if(row == 0 && column == 0) {
                    return new StaticTextField("Root");
                }
                else if(c == null) {
                    return new EmptyPanel();
                }
                else return c;
        }
    }//end of CellRenderer

    /**
     * This method sets up the property sheet for the given element.
     */
    synchronized void setTarget(Simulation.Element e) {
        updateThread.clear();
        updateThread = new UpdateThread(this);
	    removeAll();
	    depth = 0;
///        groups = new java.util.Vector();
	    // We make the panel invisivle during the reconfiguration
	    // to try to reduce screen flicker.

	    // As a workaround for #4056424, we avoid making the panel
	    // invisible first time though, during startup.
	    if (target != null)
	        setVisible(false);

	    target = e;

        //set up the table
        PropertyNode rootNode = new PropertyNode(target);
        if(rootNode == null) return;
        addChildren(rootNode);
        model = new PropertyModel(rootNode);
        treeTable = new JTreeTable(model);
        treeTable.getTree().putClientProperty("JTree.lineStyle","Angled");

        //set icons so that tree doesn't display standard folder/file icons
        ((DefaultTreeCellRenderer)treeTable.getTree().getCellRenderer()).setLeafIcon(PropertySheet.LEAF_ICON);
        ((DefaultTreeCellRenderer)treeTable.getTree().getCellRenderer()).setClosedIcon(PropertySheet.CLOSED_ICON);
        ((DefaultTreeCellRenderer)treeTable.getTree().getCellRenderer()).setOpenIcon(PropertySheet.OPEN_ICON);

        treeTable.setDefaultRenderer(Object.class, new CellRenderer());
        treeTable.setDefaultEditor(Object.class, new CellEditor());
        
        //listener that fills in the properties of an object if it is expanded by the user
        treeTable.getTree().addTreeExpansionListener(new TreeExpansionListener() {
            public void treeExpanded(TreeExpansionEvent evt) {
                PropertyNode node = (PropertyNode)evt.getPath().getLastPathComponent();
                node.removeAllChildren();
                addChildren(node);
            }
            public void treeCollapsed(TreeExpansionEvent evt) {}
        });
    
        sp = new JScrollPane(treeTable); // Put tree in a scrollable pane
        
         //Creates JComboBox that holds all currently added simulation elements.  If a different element is
         //selected from the box, the properties of the newly selected element are displayed in the prop. sheet
	    JComboBox elementCombo = new JComboBox(Simulation.instance.allElements().toArray());
	    elementCombo.setSelectedItem(target);
	    elementCombo.addItemListener(new java.awt.event.ItemListener() {
	        public void itemStateChanged(java.awt.event.ItemEvent ie){
	            Simulation.Element item = (Simulation.Element)((JComboBox)ie.getSource()).getSelectedItem();
	            Etomica.propertySheet().setTarget(item);//new Wrapper(item, item.getName(), "etomica" + item.getName()));
	        }
	    });
 
        // GridBagLayout is created and set as the current layout of the property sheet.
        GridBagConstraints gbc = new GridBagConstraints();
        GridBagLayout gbl = new GridBagLayout();
        setLayout(gbl);
        
        /*
         * Add elementCombo at (0,0) position (gridx = gridy = 0).  It can grow in x-dir but not y-dir 
         * (weightx = 1, weighty = 0).  It is allocated one gridblock in height and infinite gridblocks in 
         * length (gridheight = 1, gridwidth = REMAINDER).  It expands to fill all of its allocated
         * gridblocks in x-dir (gbc.fill = GridBagConstraints.HORIZONTAL)
         */
        gbc.gridx = gbc.gridy = 0;
        gbc.weightx = 1;
        gbc.weighty = 0;
        gbc.gridheight = 1;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbl.setConstraints(elementCombo, gbc);
        add(elementCombo);
        // end of JComboBox, elementCombo, addition
        
        /*
         * Add sp at (0,1) position (gridx = 0, gridy = 1).  It can grow in both the x-dir and the y-dir 
         * (weightx = 1, weighty = 1).  It is allocated one gridblock in height and infinite gridblocks in 
         * length (gridheight = 1, gridwidth = REMAINDER).  It expands to fill all of its allocated
         * gridblocks in x-dir and y-dir (gbc.fill = GridBagConstraints.BOTH)
         */
        gbc.gridy = 1;
        gbc.weighty = 1;
        gbc.fill = GridBagConstraints.BOTH;
        gbl.setConstraints(sp, gbc);
        add(sp);
        // end of scrollpane, sp, addition.
	    
	    frame.getContentPane().removeAll();
	    frame.getContentPane().add(this);
	    doLayout(true);

	    processEvents = true;

	    Insets ins = frame.getInsets();
	    int frameWidth = getSize().width + ins.left + ins.right + 20;
	    int frameHeight = getSize().height + ins.top + ins.bottom + 20;
        
        frameHeight += 25;//account for height of elementCombo
        frameWidth += 200;
	    frame.setSize(frameWidth,frameHeight);
        setVisible(true);
    }//end of setTarget
    
    /**
     * Fills in the properties of the given node object so they are displayed when the 
     * node is expanded.
     */
    private void addChildren(PropertyNode parentNode) {
        
        Object object = parentNode.object();
        if(object == null) return;
        
        if(object.getClass().isArray()) 
            addArrayChildren(parentNode);
        else 
            addObjectChildren(parentNode);
            
        //notify model of the changes below this node
	    model.nodeStructureChanged(parentNode);
	    
    }
    
    //called by addChildren
    //sets up children of an array property
    private void addArrayChildren(PropertyNode parentNode) {
        
        Object[] objects = (Object[])parentNode.object();
                
        for(int i=0; i<objects.length; i++) {
            Object object = objects[i];
            Object value = object;
            MyLabel newLabel = new MyLabel("["+i+"]", Label.LEFT);
            Component view = new EmptyPanel(value.toString());
            PropertyNode child = new PropertyNode(value, newLabel, view, null, null, null);
            
//See if object can have child objects in tree (cannot if it is primitive)
            if(!(value == null || 
                value instanceof Number || 
                value instanceof Boolean ||
                value instanceof Character ||
                value instanceof String ||
                value instanceof Color ||
                value instanceof etomica.Constants.TypedConstant ||
                value instanceof java.awt.Font)) {/*add dummy child*/
                child.add(new PropertyNode(null,new JLabel(),new EmptyPanel(), new EmptyPanel(), null, null));
            }
            
            parentNode.add(child);
            System.out.println(((etomica.MCMove)object).getName());
        }
        

        //Introspection to get array of all properties
//        PropertyDescriptor properties = new PropertyDescriptor;
    }
    
    //called by addChildren
    //sets up children of a non-array property
    private void addObjectChildren(PropertyNode parentNode) {
        
        Object object = parentNode.object();

        //Introspection to get array of all properties
        PropertyDescriptor[] properties = null;
        BeanInfo bi = null;
        try {
	        bi = Introspector.getBeanInfo(object.getClass());
	        properties = bi.getPropertyDescriptors();
	    } 
	    catch (IntrospectionException ex) {
	        error("PropertySheet: Couldn't introspect", ex);
	        return;
	    }

        //Loop through all properties and determine current value and find appropriate editor
        DimensionedDoubleEditor averageEditor = null;
        DimensionedDoubleEditor errorEditor = null;
        DimensionedDoubleEditor currentEditor = null;
        PropertyText averageView = null;
        PropertyText errorView = null;
        PropertyText currentView = null;
	    for (int i = 0; i < properties.length; i++) {
	        PropertyNode node = processProperty(parentNode, properties[i], bi);
	        if(node == null) continue;
	        String name = properties[i].getDisplayName();  //Localized display name
	        //if meter, set up update thread that refreshes meter's measured values
    	    if(object instanceof Meter) {
		        if(name.equals("average")) {
		            averageEditor = (DimensionedDoubleEditor)node.editor();
		            averageView = (PropertyText)node.view();
		        }
		        else if(name.equals("error")) {
		            errorEditor = (DimensionedDoubleEditor)node.editor();
		            errorView = (PropertyText)node.view();
		        }
		        else if(name.equals("mostRecent")) {
		            currentEditor = (DimensionedDoubleEditor)node.editor();
		            currentView = (PropertyText)node.view();
		        }
	        }//end if
	    }//end of loop over properties

	    //if meter, set to display/update its values on table
	    if(object instanceof Meter) updateThread.add((Meter)object,currentEditor,averageEditor,errorEditor,
	                                        currentView,  averageView,  errorView  );
	                                        
    }//end of addObjectChildren method
    
    private PropertyNode processProperty(PropertyNode parentNode, PropertyDescriptor property,
                                            BeanInfo bi) {
        
	        // Don't display hidden or expert properties.
	        if (property.isHidden() || property.isExpert())
		        return null;

	        Object value = null;
	        Component view = null;
	        Component unitView = null;
	        JLabel label = null;
	        PropertyEditor editor = null;

	        String name = property.getDisplayName();  //Localized display name 
	        if(name.equals("dimension")) return null;         //skip getDimension()
	        Class type = property.getPropertyType();  //Type (class) of this property
	        Method getter = property.getReadMethod(); //method used to read value of property in this object
	        Method setter = property.getWriteMethod();//method used to set value of property
	        // Only display read/write properties.
	        if (getter == null || (setter == null && !(parentNode.object() instanceof Meter)) || type == Class.class) {
		        return null;
	        }
	        try {
	            //read the current value of the property
		        Object args[] = { };
		        try {value = getter.invoke(parentNode.object(), args);}
		        catch(NullPointerException ex) {value = null;}

                //find and instantiate the editor used to modify value of the property
	            if(property.isConstrained())
	                editor = new ConstrainedPropertyEditor();
	            //if property is a TypedConstant
	            else if(etomica.Constants.TypedConstant.class.isAssignableFrom(type) && value != null) {
	                editor = new TypedConstantEditor();
	            }
	            else if(etomica.units.Unit.class.isAssignableFrom(type)) {
	                editor = new etomica.UnitEditor((Unit)value);
	            }
	            else {
        	        //property is a dimensioned number
        	        if(type == Double.TYPE) {
        	            //try to get dimension from get(property)Dimension() method
        	            etomica.units.Dimension dimension = etomica.units.Dimension.introspect(parentNode.object(),name,bi);
        	            //try to get dimension from getDimension() method
                        if(dimension == null) dimension = etomica.units.Dimension.introspect(parentNode.object(),"",bi);
        	            if(dimension != null) {
        	                editor = new DimensionedDoubleEditor(dimension);
        	            }
        	        }
        	        //property is not a dimensioned number; see if its editor was set explicitly
        	        if(editor == null) { 
	                    Class pec = property.getPropertyEditorClass();
		                if (pec != null) {
		                    try {
			                    editor = (PropertyEditor)pec.newInstance();
		                    } 
		                    catch (Exception ex) {}
		                }
		            }
		            //property is not a dimensioned number and was not set explicitly
		            //have editor manager look for an appropriate editor
		            if (editor == null)
		                editor = PropertyEditorManager.findEditor(type);
		        }//done with trying to get an editor for the property
		        
	            // If we can't edit this component, skip it.
	            if (editor == null) {
		            // If it's a user-defined property we give a warning.
		            String getterClass = property.getReadMethod().getDeclaringClass().getName();
		            if (getterClass.indexOf("java.") != 0) {
		                System.err.println("Warning: Can't find public property editor for property \"" + name + "\".  Skipping.");
		            }
		            return null;
	            }

                //set the editor to the current value of the property
	            try {
	                editor.setValue(value);
	            } catch(NullPointerException e) {}
	            
	            //add listener that causes the wasModified method to be 
	            //invoked when editor fires property-change event
	            editor.addPropertyChangeListener(adaptor);

		    // Now figure out how to display it...
		        if (editor.isPaintable() && editor.supportsCustomEditor())
        		    view = new PropertyCanvas(frame, editor);
        		else if (editor instanceof etomica.UnitEditor) 
        		    view = ((etomica.UnitEditor)editor).unitSelector();
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
		            return null;
		        }
		        if(editor instanceof DimensionedDoubleEditor) {
		            unitView = ((DimensionedDoubleEditor)editor).unitSelector();
		            if(parentNode.object() instanceof etomica.Meter)
		                ((PropertyText)view).setEditable(false);
		        }
		        else unitView = new EmptyPanel();
	        } //end of try
	        catch (InvocationTargetException ex) {
		        System.err.println("Skipping property " + name + " ; exception on target: " + ex.getTargetException());
		        ex.getTargetException().printStackTrace();
		        return null;
	        } 
	        catch (Exception ex) {
		        System.err.println("Skipping property " + name + " ; exception: " + ex);
		        ex.printStackTrace();
		        return null;
	        }

            MyLabel newLabel = new MyLabel(name, Label.LEFT);
	        
            PropertyNode child = new PropertyNode(value, newLabel, view, unitView, editor, property);
            
//See if object can have child objects in tree (cannot if it is primitive)
            if(!(value == null || 
                value instanceof Number || 
                value instanceof Boolean ||
                value instanceof Character ||
                value instanceof String ||
                value instanceof Color ||
                value instanceof etomica.Constants.TypedConstant ||
                value instanceof java.awt.Font)) {/*add dummy child*/
                child.add(new PropertyNode(null,new JLabel(),new EmptyPanel(), new EmptyPanel(), null, null));
            }
            
            parentNode.add(child);
            
            return child;
    }//end of processProperty
        
    private void doLayout(boolean doSetSize) {
        JTree tree = treeTable.getTree();
        PropertyNode root = (PropertyNode)tree.getModel().getRoot();
        if(root==null || root.getChildCount()==0) return;
      //  if(groups.isEmpty())  return;

	    // First figure out the size of the columns.
	    int labelWidth = 92;
	    int viewWidth = 120;

        for(java.util.Enumeration enum = root.breadthFirstEnumeration(); enum.hasMoreElements();) {
            PropertyNode node = (PropertyNode)enum.nextElement();
            if(node == root || node.label() == null || node.view() == null) continue;

	        int w = node.label().getPreferredSize().width;
	        if (w > labelWidth) {
		        labelWidth = w;
	        }
	        w = node.view().getPreferredSize().width;
	        if (w > viewWidth) {
		        viewWidth = w;
	        }
	    }//end for
	    int width = 3*hPad + labelWidth + viewWidth;

	    // Now position all the components.
	    int y = 10;
        for(java.util.Enumeration enum = root.breadthFirstEnumeration(); enum.hasMoreElements();) {
            PropertyNode node = (PropertyNode)enum.nextElement();
            if(node == root || node.label() == null || node.view() == null) continue;

	        node.label().setBounds(hPad, y, labelWidth, 20);
	        int h = 20;
	        node.view().setBounds(labelWidth + 2*hPad, y, viewWidth, h);
	        y += (h + vPad);
	    }
	    y += vPad;

	    if(doSetSize) setSize(width, y);
    }
    
/*    private void doLayout(boolean doSetSize) {
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
	        int h = 20;
	        g.views[i].setBounds(labelWidth + 2*hPad, y, viewWidth, h);
	        y += (h + vPad);
	    }
	    y += vPad;

	    if (doSetSize)
	        setSize(width, y);
    }
*/
    /**
     * When property event is fired, this method transmits the change to the actual
     * instance of the object being edited.  Also checks for change in value of any
     * other properties and repaints the sheet accordingly.
     */
    synchronized void wasModified(PropertyChangeEvent evt) {
        
        if (!processEvents) return;

        JTree tree = treeTable.getTree();
        PropertyNode root = (PropertyNode)tree.getModel().getRoot();
	    if (evt.getSource() instanceof PropertyEditor) {
	        PropertyEditor editor = (PropertyEditor) evt.getSource();
            //loop over all objects and properties to find the editor that made the change
            //***could replace this process with a hash table***
            for(java.util.Enumeration enum = root.breadthFirstEnumeration(); enum.hasMoreElements();) {
                PropertyNode node = (PropertyNode)enum.nextElement();
                
                if(node==root || node.editor()!=editor) continue;
                
                //found editor
		        PropertyDescriptor property = node.descriptor();
		        Object value = editor.getValue();
		        
		        //implement change in value stored in node
		        node.setObject(value);
		        
		        Method setter = property.getWriteMethod();
		        if(setter == null) continue;
		        
		        //write java code that implements change for writing to file, if requested later
		        JavaWriter javaWriter = (JavaWriter)Etomica.javaWriters.get(target.parentSimulation());
                javaWriter.propertyChange(target, setter, editor.getJavaInitializationString());
		        
		        //implement the change in the object being edited
		        try {
		            Object args[] = { value };
		            args[0] = value;
		            setter.invoke(((PropertyNode)node.getParent()).object(), args);
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
	  //              node.view().repaint();
	            }
	            if(node.getChildCount() > 0) {
                    node.removeAllChildren();
	                addChildren(node);
	            }
		        break;//found editor; stop looking
	        }//end of for(Enumeration loop
	    }

	    // Now re-read all the properties and update the editors
	    // for any other properties that have changed.
	    Object args[] = { };
        for(java.util.Enumeration enum = root.breadthFirstEnumeration(); enum.hasMoreElements();) {
            PropertyNode node = (PropertyNode)enum.nextElement();
            if(node == root) continue;
            Object target = ((PropertyNode)node.getParent()).object();
	        Object o;
	        Method setter = null;
	        try {
	            Method getter = node.descriptor().getReadMethod();
	            setter = node.descriptor().getWriteMethod();
	            o = getter.invoke(target, args);
	        } 
	        catch (Exception ex) { o = null; }
	        if (o == node.object() || (o != null && o.equals(node.object())) ||  setter == null 
	            || (o.getClass().isArray() && node.object().getClass().isArray() && java.util.Arrays.equals((Object[])o,(Object[])node.object()))) {
	            // The property is equal to its old value.
    		    continue;
	        }
	        node.setObject(o);
	        // Make sure we have an editor for this property...
	        if (node.editor() == null)
		        continue;
    	        
	        // The property has changed!  Update the editor.
	        node.editor().setValue(o);
	        if (node.view() != null)
		        node.view().repaint();

	        // Make sure the target bean gets repainted.
	        if (Beans.isInstanceOf(target, Component.class)) {
	            ((Component)(Beans.getInstanceOf(target, Component.class))).repaint();
	        }
	    }//end of for(java.util.Enumeration) loop
	        
	    //refresh the display of the simulation
        ((SimulationFrame)Etomica.simulationFrames.get(Simulation.instance)).repaint();
        
        //refresh the display of the property sheet
        //this handles repainting of values when units change, and unitlist choices when
        //prefix changes.
        //unresolved problem: units in unitlist box that spill outside the propertysheet
        //                    frame are not repainted with new prefix
        frame.repaint();
	    
    }//end of wasModified method

   
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
    
/*    private static class PropertyGroup {
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
*/    
    private static class EmptyPanel extends JTextField {
        EmptyPanel() {
            super("");
            setBorder(PropertySheet.EMPTY_BORDER);
        }
        EmptyPanel(String text) {
            super(text);
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
    private Simulation.Element target;
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
