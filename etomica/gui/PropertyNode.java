package etomica.gui;
import etomica.units.Unit;

import java.awt.Component;
import javax.swing.JLabel;

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
    
    public PropertyNode() {
        super();
    }
    public PropertyNode(JLabel newLabel, Component newView, Component newUnitView) {
        super(newLabel);
        label = newLabel;
        view = newView;
        unitView = newUnitView;
    }
    
    public JLabel label() {return label;}
    public Component view() {return view;}
    public Component unitView() {return unitView;}
    
}