package etomica.lattice;
import etomica.*;

/**
 * Iterates over the neighbors of a particular site, as specified by 
 * the site's neighborManager.
 */
public class SiteIteratorNeighbor implements AtomIterator {
    
    private NeighborManager neighborManager;
    private final AtomIteratorList iterator = new AtomIteratorList();
    private boolean upListNow, doGoDown;
    private Atom next;
    
    public boolean hasNext() {return next != null;}
    public Atom reset() {
        throw new RuntimeException("method SiteIteratorNeighbor.reset() not yet implemented");
    }
    public void setBasis(Atom atom) {
        throw new RuntimeException("method SiteIteratorNeighbor.setBasis(Atom) not yet implemented");
    }
    public Atom getBasis() {
        throw new RuntimeException("method SiteIteratorNeighbor.getBasis() not yet implemented");
    }
    public void setAsNeighbor(boolean b) {
        throw new RuntimeException("method SiteIteratorNeighbor.setAsNeighbor not implemented");
    }
    public Atom reset(IteratorDirective d) {
        throw new RuntimeException("method SiteIteratorNeighbor.reset(IteratorDirective) not yet implemented");
    }
    public boolean contains(Atom a) {
        throw new RuntimeException("method SiteIteratorNeighbor.contains(Atom) not yet implemented");
    }
        
    public void reset(Atom site, IteratorDirective.Direction direction) {
        reset(((Site)site).neighborManager(), direction);
    }
    public void reset(NeighborManager manager, IteratorDirective.Direction direction) {
        neighborManager = manager;
        if(direction == null) direction = IteratorDirective.NEITHER;
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
        iterator.reset((AtomLinker)null);
        next = null;
        if(neighborManager == null) return;
        if(upListNow) iterator.setBasis(neighborManager.upNeighbors());
        if(iterator.hasNext()) {
            next = iterator.next();
        } else if(doGoDown) {
            iterator.setBasis(neighborManager.downNeighbors());
            if(iterator.hasNext()) next = iterator.next();
            doGoDown = false;
        }
    }
    public Atom first() {
        throw new RuntimeException("method first() not implemented in SiteIteratorNeighbor");
    }
    public Atom next() {
        Atom nextAtom = next;
        if(iterator.hasNext()) {
            next = iterator.next();
        } else if(doGoDown) {
            iterator.setBasis(neighborManager.downNeighbors());
            
            if(iterator.hasNext()) next = iterator.next();
            else next = null;
            
            doGoDown = false;
        } else next = null;
        
        return nextAtom;
    }
    public void allAtoms(AtomAction act) {
        if(neighborManager == null) return;
        if(upListNow) {
            iterator.setBasis(neighborManager.upNeighbors());
            iterator.allAtoms(act);
        }
        if(doGoDown) {
            iterator.setBasis(neighborManager.downNeighbors());
            iterator.allAtoms(act);
        }
    }
    
    public int size() {return neighborManager.neighborCount();}
}