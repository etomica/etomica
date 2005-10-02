package etomica.util;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.atom.SpeciesRoot;

/**
 * Extends the java ObjectInputStream to add fields and methods to enable complete
 * reconstruction of simulation from a serialized form.
 */

/*
 * Created on Jul 6, 2005
 */

public class EtomicaObjectInputStream extends ObjectInputStream {

    public SpeciesRoot speciesRoot;
    public LinkedList dirtyObjects = new LinkedList();
    public HashMap objectData = new HashMap();
    
    /**
     * @param out
     * @throws java.io.IOException
     */
    public EtomicaObjectInputStream(InputStream out) throws IOException {
        super(out);
    }

    /**
     * @throws java.io.IOException
     * @throws java.lang.SecurityException
     */
    public EtomicaObjectInputStream() throws IOException, SecurityException {
        super();
    }

    public void finalizeRead() {
        Iterator listIterator = dirtyObjects.iterator();
        while (listIterator.hasNext()) {
            DirtyObject obj = (DirtyObject)listIterator.next();
            obj.rebuild(objectData.get(obj));
        }
        dirtyObjects.clear();
        objectData.clear();
    }
}
