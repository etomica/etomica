package simulate.gui;

import simulate.gui.treetable.*;

public class PropertyModel extends AbstractTreeTableModel {
    
    PropertyModel(PropertyNode root) {
        super(root);
    }
        
    public Object getChild(Object parent, int index) {
        return ((PropertyNode)parent).getChildAt(index);
    }
        
    public int getChildCount(Object parent) {
        return ((PropertyNode)parent).getChildCount();
    }
        
    public int getColumnCount() {
        return 3;
    }
        
    public Class getColumnClass(int column) {
        switch(column) {
            case 0: return TreeTableModel.class;
            case 1: return Object.class;
            case 2: return Object.class;
            default: return null;
        }
    }

    public String getColumnName(int column) {
        switch(column) {
            case 0: return "Property";
            case 1: return "Value";
            case 2: return "Units";
            default: return null;
        }
    }
        
    public Object getValueAt(Object node, int column) {
        PropertyNode propertyNode = (PropertyNode)node;
        switch(column) {
            case 0: return propertyNode.label();
            case 1: return propertyNode.view();
            case 2: return propertyNode.unitView();
            default: return null;
        }
    }
    
   /** Override superclass so that by default, the column without the Tree in it the only editable one. 
    *  Making this column editable causes the JTable to forward mouse 
    *  and keyboard events in the Tree column to the underlying JTree. 
    */ 
    public boolean isCellEditable(Object node, int column) { 
         return true; 
    }

}
