/*
 * Created on Aug 23, 2004 
 */
package etomica.utility;

import etomica.Debug;

/**
 * Linker used in the binary tree TreeList.  TreeLinkers point to the
 * TreeLinker above, to the left and right of them within the tree.
 * 
 * @author andrew
 */
public class TreeLinker {
    public final Object obj;
    public TreeLinker previous, nextLeft, nextRight;
    public double sortKey;

    /**
     * Constructor throws exception if given atom is null.  Only
     * AtomLink.Tab instances can have a null atom field.
     */
    public TreeLinker(Object o) {
        obj = o;
        sortKey = 0.0;
        if (Debug.ON) {
            nextLeft = previous = nextRight = null;
        }
    }
    
    public String toString() {
        return "TreeLinker ("+obj.toString()+" "+sortKey+")";
    }

    /**
     * Disconnects this linker from the linkers above and below it, and
     * puts them in sequence, repairing the hole.  There is no indication 
     * in the linker itself that it is no longer part of the tree (unless 
     * Debug is on).
     */
    public void remove() {
        TreeLinker newNext;
        if (nextRight == null) {
            newNext = nextLeft;
        }
        else if (nextLeft == null) {
            newNext = nextRight;
        }
        else {
            if (nextRight.nextLeft == null) {
                newNext = nextRight;
            }
            else {
                // this section of the algorithm is not meant to be understood
                // it is meant to work.  So save the Tylenol and believe.
                newNext = nextRight.nextLeft;
                while (newNext.nextLeft != null) {
                    newNext = newNext.nextLeft;
                }

                if (newNext.nextRight != null) {
                    newNext.nextRight.previous = newNext.previous;
                }
                newNext.previous.nextLeft = newNext.nextRight;
                nextRight.previous = newNext;
                newNext.nextRight = nextRight;
            }
            nextLeft.previous = newNext;
            newNext.nextLeft = nextLeft;
        }

        if (newNext != null) {
            newNext.previous = previous;
        }

        if (previous.nextLeft == this) {
            previous.nextLeft = newNext;
        }
        else {
            previous.nextRight = newNext;
        }
        if (Debug.ON) {
            previous = nextRight = nextLeft = null;
        }
    }
}
