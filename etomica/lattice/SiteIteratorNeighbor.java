package etomica.lattice;
import etomica.IteratorDirective;

/**
 * Iterates over the neighbors of a particular site, as specified by 
 * the site's neighborManager.
 */
public class SiteIteratorNeighbor implements SiteIterator {
    
    private NeighborManager neighborManager;
    private final SiteIteratorList iterator = new SiteIteratorList();
    private boolean upListNow, doGoDown;
    private Site next;
    
    public boolean hasNext() {return next != null;}
    public void reset() {
    }
    public void reset(Site site, IteratorDirective.Direction direction) {
        reset(site.neighborManager(), direction);
    }
    public void reset(NeighborManager manager, IteratorDirective.Direction direction) {
        neighborManager = manager;
        if(direction == null) direction = IteratorDirective.NEITHER;
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
        iterator.reset((SiteLinker)null);
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
    public Site first() {
        throw new RuntimeException("method first() not implemented in SiteIteratorNeighbor");
    }
    public Site next() {
        Site nextSite = next;
        if(iterator.hasNext()) {
            next = iterator.next();
        } else if(doGoDown) {
            iterator.setBasis(neighborManager.downNeighbors());
            
            if(iterator.hasNext()) next = iterator.next();
            else next = null;
            
            doGoDown = false;
        } else next = null;
        
        return nextSite;
    }
    public void allSites(SiteAction act) {
        if(neighborManager == null) return;
        if(upListNow) {
            iterator.setBasis(neighborManager.upNeighbors());
            iterator.allSites(act);
        }
        if(doGoDown) {
            iterator.setBasis(neighborManager.downNeighbors());
            iterator.allSites(act);
        }
    }
    
    public int size() {return neighborManager.neighborCount();}
}