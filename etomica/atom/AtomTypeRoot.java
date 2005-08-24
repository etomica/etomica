package etomica.atom;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jul 11, 2005 by kofke
 */
public class AtomTypeRoot extends AtomTypeGroup {

    /**
     * Used only to create root type
     */
    AtomTypeRoot(AtomIndexManager indexManager) {
        super(indexManager);
    }

    int requestIndex() {
        return ++descendantCount;
    }
    
    private int descendantCount = 0;
    
}
