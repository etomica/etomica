/*
 * Created on Aug 23, 2004
 */
package etomica.utility;

import etomica.Debug;

/**
 * TreeList is a generic binary tree structure made up of TreeLinkers.
 * Objects are added with a sortKey, which determines how it is added
 * to the tree (left for objects with lower sort keys, right for objects
 * with higher sort keys).  The lowest event is always at the bottom on
 * the left.  Objects can be added and removed.
 *
 * @author andrew
 */

/* The treelist is powerful, but is not very hard to screw up.  It can
 * get circular links among other things.  Do not
 * add an object (linker) which is already in the tree. Do not remove
 * an object that is not in the tree.  Do not taunt the tree.  Turning
 * on Debug will attempt to catch problems.
 */
public class TreeList {

    private TreeLinker head;

    public TreeList() {
        head = new TreeLinker(null);
    }

    public void add(TreeLinker newNode) {
        if (Debug.ON && (newNode.nextLeft != null || newNode.nextRight != null || newNode.nextLeft != null)) {
            throw new RuntimeException("attempting to add a node to the tree which is already part of a tree");
        }
        TreeLinker parentNode = head;
        if (head.nextRight == null) {
            head.nextRight = newNode;
            newNode.previous = head;
        }
        else {
            while (true) {
                if (newNode.sortKey > parentNode.sortKey) {
                    if (parentNode.nextRight == null) {
                        parentNode.nextRight = newNode;
                        newNode.previous = parentNode;
                        break;
                    }
                    parentNode = parentNode.nextRight;
                    continue;
                }
                if (parentNode.nextLeft == null) {
                    parentNode.nextLeft = newNode;
                    newNode.previous = parentNode;
                    break;
                }
                parentNode = parentNode.nextLeft;
            }
        }
        newNode.nextRight = newNode.nextLeft = null;
    }

    static public void remove(TreeLinker oldNode) {
        oldNode.remove();
    }

    /**
     * returns the object from the binary tree associated
     * with the lowest sort-key, or null if the tree is empty
     * @return first object in the tree
     */
    public Object firstNode() {
        TreeLinker first = head.nextLeft;
        if (first==null) {
            first = head.nextRight;
        }
        while (first.nextLeft != null) {
            first = first.nextLeft;
        }
        if (first != null) {
            return first.obj;
        }
        return null;
    }

    /**
     * clear the tree by disconnecting everything from the head.  Individual
     * linkers will still be linked to each other (unless Debug is ON).
     */
    public void reset() {
        if (Debug.ON) {
            TreeLinker node;
            while ((node = (TreeLinker)firstNode()) != null) {
                node.remove();
            }
        }
        head.nextLeft = null;
        head.nextRight = null;
    }
}
