package etomica.gui;

import etomica.gui.treetable.*;
import javax.swing.tree.TreeNode;

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

    
    //the following methods are copied from DefaultTreeModel
    
    /**
      * Invoke this method if you've totally changed the children of
      * node and its childrens children...  This will post a
      * treeStructureChanged event.
      */
    public void nodeStructureChanged(TreeNode node) {
        if(node != null) {
           fireTreeStructureChanged(this, getPathToRoot(node), null, null);
        }
    }
    
    /**
     * Builds the parents of node up to and including the root node,
     * where the original node is the last element in the returned array.
     * The length of the returned array gives the node's depth in the
     * tree.
     * 
     * @param aNode the TreeNode to get the path for
     * @param an array of TreeNodes giving the path from the root to the
     *        specified node. 
     */
    public TreeNode[] getPathToRoot(TreeNode aNode) {
        return getPathToRoot(aNode, 0);
    }
    
    /**
     * Builds the parents of node up to and including the root node,
     * where the original node is the last element in the returned array.
     * The length of the returned array gives the node's depth in the
     * tree.
     * 
     * @param aNode  the TreeNode to get the path for
     * @param depth  an int giving the number of steps already taken towards
     *        the root (on recursive calls), used to size the returned array
     * @return an array of TreeNodes giving the path from the root to the
     *         specified node 
     */
    protected TreeNode[] getPathToRoot(TreeNode aNode, int depth) {
        TreeNode[]              retNodes;
	// This method recurses, traversing towards the root in order
	// size the array. On the way back, it fills in the nodes,
	// starting from the root and working back to the original node.

        /* Check for null, in case someone passed in a null node, or
           they passed in an element that isn't rooted at root. */
        if(aNode == null) {
            if(depth == 0)
                return null;
            else
                retNodes = new TreeNode[depth];
        }
        else {
            depth++;
            if(aNode == root)
                retNodes = new TreeNode[depth];
            else
                retNodes = getPathToRoot(aNode.getParent(), depth);
            retNodes[retNodes.length - depth] = aNode;
        }
        return retNodes;
    }
}
