package etomica.gui;
import etomica.units.Unit;

import java.awt.Component;
import javax.swing.JLabel;
import java.beans.PropertyEditor;
import java.beans.PropertyDescriptor;

/**
 * Node on a property tree for an object.  This node type extends the DefaultMutableTreeNode
 * of tree to permit it to hold a label and a Component.  The component should be the
 * graphical interface of a property editor that implements changes to the value of the property.
 *
 * @see PropertySheet
 */

public class PropertyNode extends javax.swing.tree.DefaultMutableTreeNode {
    
    private JLabel label;
    private Component view;
    private Component unitView;
    private Unit unit;
    private Object object;
    private PropertyEditor editor;
    private PropertyDescriptor descriptor;
    
    /**
     * Constructor for the root node, corresponding to the object with the properties being edited.
     */
    public PropertyNode(Object newObject) {
        super(newObject!=null ? newObject.toString() : "No selection");
        object = newObject;
    }
    /**
     * Constructor for a property node of the object being edited.
     */
    public PropertyNode(Object newObject, 
                        JLabel newLabel, Component newView, Component newUnitView,
                        PropertyEditor ed, PropertyDescriptor desc) {
        super(newLabel);
        object = newObject;
        label = newLabel;
        view = newView;
        unitView = newUnitView;
        editor = ed;
        descriptor = desc;
    }
   
    public Object object() {return object;}
    public void setObject(Object obj) {object = obj;}
    public JLabel label() {return label;}
    public Component view() {return view;}
    public Component unitView() {return unitView;}
    public PropertyEditor editor() {return editor;}
    public PropertyDescriptor descriptor() {return descriptor;}
    
}