package etomica.util;


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
public class TreeList implements java.io.Serializable {

    private static final long serialVersionUID = 1L;
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
            parentNode = parentNode.nextRight;
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
    public Object firstElement() {
        TreeLinker first = head.nextLeft;
        if (first==null) {
            first = head.nextRight;
            if (first == null) {
                return null;
            }
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
     * returns the object from the binary tree associated
     * with the lowest sort-key, or null if the tree is empty
     * @return first object in the tree
     */
    private TreeLinker firstNode() {
        TreeLinker first = head.nextLeft;
        if (first==null) {
            first = head.nextRight;
            if (first == null) {
                return null;
            }
        }
        while (first.nextLeft != null) {
            first = first.nextLeft;
        }
        if (first != null) {
            return first;
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
            while ((node = firstNode()) != null) {
                node.remove();
            }
        }
        head.nextLeft = null;
        head.nextRight = null;
    }
    
    /**
     * checks the tree for consistency of connections and ordering by sortKey.
     * circular (but consistent) links will cause this method to loop forever.
     */
    public void check() {
        TreeLinker node = head;
        if (head.nextLeft != null || head.previous != null) {
            System.out.println("head: p "+head.previous+" l "+head.nextLeft+" r "+head.nextRight);
            throw new IllegalStateException("head's previous and nextRight must be null");
        }
        node = head.nextRight;
        if (node == null) {
            // tree is empty
            return;
        }
        if (node.previous != head) {
            System.out.println("head.r "+head.nextRight+" head.r.p "+head.nextRight.previous);
            throw new IllegalStateException("head.r.p is not head");
        }
        while (true) {
            // check this node
            if (node.previous.nextRight != node && node.previous.nextLeft != node) {
                System.out.println("node "+node+" node.p.r "+node.previous.nextRight+" node.p.l "+node.previous.nextLeft);
                throw new IllegalStateException("connections inconsistent");
            }
            if (node.previous != head) {
                if (node.previous.nextRight == node && node.previous.sortKey > node.sortKey) {
                    System.out.println("node "+node+" node.sortKey "+node.sortKey+" node.p "+node.previous+" node.p.sortKey "+node.previous.sortKey);
                    throw new IllegalStateException("node is right from parent and should have a higher sortKey");
                }
                else if (node.previous.nextLeft == node && node.previous.sortKey < node.sortKey) {
                    System.out.println("node "+node+" node.sortKey "+node.sortKey+" node.p "+node.previous+" node.p.sortKey "+node.previous.sortKey);
                    throw new IllegalStateException("node is left from parent and should have a lower sortKey");
                }
            }
            if (node.nextLeft != null) {
                if (node.nextLeft.previous != node) {
                    System.out.println("node "+node+" node.l "+node.nextLeft+" node.l.p "+node.nextRight.previous);
                    throw new IllegalStateException("connections inconsistent");
                }
                if (node.sortKey < node.nextLeft.sortKey) {
                    System.out.println("node "+node+" node.l "+node.nextLeft);
                    throw new IllegalStateException("node's nextLeft should have a lower sortKey");
                }
            }
            if (node.nextRight != null) {
                if (node.nextRight.previous != node) {
                    System.out.println("node "+node+" node.r "+node.nextRight+" node.l.p "+node.nextRight.previous);
                    throw new IllegalStateException("connections inconsistent");
                }
                if (node.sortKey > node.nextRight.sortKey) {
                    System.out.println("node "+node+"\nnode.r "+node.nextRight);
                    throw new IllegalStateException("node's nextRight should have a higher sortKey");
                }
            }
            if (node.nextLeft != null) {
                // go down left
                node = node.nextLeft;
//                System.out.println("down left 0");
                continue;
            }
            if (node.nextRight != null) {
                // go down right
                node = node.nextRight;
//                System.out.println("down right 0");
                continue;
            }
            if (node.previous.nextLeft == node) {
                // go back up right until we can go down right
                while (node.nextRight == null && node != head && node.previous.nextLeft == node) {
                    node = node.previous;
//                    System.out.println("up right 0");
                }
                if (node == head) {
                    // head.nextRight.nextRight was null
                    // success
                    return;
                }
                if (node.nextRight != null) {
                    // go down right
//                    System.out.println("down right 1");
                    node = node.nextRight;
                    continue;
                }
                // can't go up right any more.  go up left instead and then up right, and then down right.
            }
            // keep going up right until we can go down right.
            while (node.nextRight == null && node != head) {
                while (node.previous.nextRight == node && node.previous != head) {
                    // up left until we can go up right
//                    System.out.println("up left 0");
                    node = node.previous;
                }
                if (node.previous == head) {
                    // success
                    return;
                }
//                System.out.println("up right 1");
                node = node.previous;
            }
            // go down right
//            System.out.println("down right 2");
            node = node.nextRight;
        }
    }
}
