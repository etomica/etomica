package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomSet;
import etomica.action.AtomsetAction;
import etomica.atom.AtomList;
import etomica.utility.Arrays;

/**
 * Iterates over all the atoms given across an array of iterators.  
 *
 * @author David Kofke
 */
 
public final class AtomIteratorListCompound implements AtomIterator {
    
    /**
     * Construct iterator to loop over all iterates obtained from each iterator 
     * in the given array.
     */
    public AtomIteratorListCompound(AtomList[] lists) {
        listSet = lists;
        iterator = new AtomIteratorListSimple();
        reset();
    }

    public AtomIteratorListCompound() {
        iterator = new AtomIteratorListSimple();
        unset();
    }
    
    /**
     * Adds the given iterator to the set of iterators collected by this 
     * compound iterator.  No check is made the iterator was not already added;
     * if an iterator is added more than once, it will be iterated for each addition,
     * as if it were a different instance each time.
     */
    public void addIterator(AtomList list) {
        if(list == null) return;
        listSet = (AtomList[])Arrays.addObject(listSet, list);
        unset();
    }
    
    public void setLists(AtomList[] lists) {
        listSet = lists;
        reset();
    }
    
    /**
     * Removes the given iterator from the set of iterators collected by
     * this compound iterator.  If iterator was not previously added (or is
     * null), no action is taken.
     */
    public void removeList(AtomList list) {
        listSet = (AtomList[])Arrays.removeObject(listSet, list);
    }
    
    public boolean hasNext() {return hasNext;}
    
    public void unset() {hasNext = false;}
    
    public int nBody() {return 1;}
    
    public int size() {
        if(listSet == null) return 0;
        int count = 0;
        for(int i=0; i<listSet.length; i++) {
            count += listSet[i].size();
        }
        return count;
    }
    
    public boolean contains(AtomSet atoms) {
        if(listSet == null) return false;
        for(int i=0; i<listSet.length; i++) {
            if(listSet[i].contains((Atom)atoms)) return true;
        }
        return false;
    }
        
    public void reset() {
        hasNext = false;
        if(!hasList()) {
            return;
        }
        index = 0;
        
        while(listSet[index].isEmpty() && index+1 < listSet.length) {
        	index++;
        }
        if (index < listSet.length) {
            iterator.setList(listSet[index]);
            iterator.reset();

            hasNext = true;
        }
    }
    
    public AtomSet peek() {
    	return hasNext ? iterator.peek() : null;
    }
    
    public AtomSet next() {
        return nextAtom();
    }//end of next

    public Atom nextAtom() {
        if(!hasNext) return null;
        Atom atom = iterator.nextAtom();
        while(!iterator.hasNext()) {
            if(++index < listSet.length) {
                iterator.setList(listSet[index]);
                iterator.reset();
            }
            else {
                hasNext = false;
                break;
            }
        }
        return atom;
    }
    
    public void allAtoms(AtomsetAction action) {
        for(int i=0; i<listSet.length; i++) {
            iterator.setList(listSet[i]);
            iterator.reset();
            iterator.allAtoms(action);
        }
    }

    private boolean hasList() {
        return (listSet.length > 0);
    }
    
    private AtomList[] listSet;
    private AtomIteratorListSimple iterator;
    private boolean hasNext;
    private int index;
    
}//end of AtomIteratorCompound