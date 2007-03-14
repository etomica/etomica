package etomica.util;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.atom.Atom;

/**
 * Extends the java ObjectInputStream to add fields and methods to enable complete
 * reconstruction of simulation from a serialized form.
 */

/*
 * Created on Jul 6, 2005
 */

public class EtomicaObjectInputStream extends ObjectInputStream {

    public LinkedList dirtyObjects = new LinkedList();
    public HashMap objectData = new HashMap();
    private static ArrayList atomList = new ArrayList();
    
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
        //sort the atom list so atoms can be found by index quickly.
        Collections.sort(atomList);
        Iterator listIterator = dirtyObjects.iterator();
        while (listIterator.hasNext()) {
            DirtyObject obj = (DirtyObject)listIterator.next();
            obj.rebuild(objectData.get(obj));
        }
        // clean up after ourselves
        dirtyObjects.clear();
        objectData.clear();
        atomList.clear();
    }
    
    public void addAtom(Atom a) {
        atomList.add(a);
    }

    /**
     * Finds the child of the given speciesMaster having the given index.
     * @throws IllegalArgumentException if no Atom is found.
     */
    public static Atom getAtomForIndex(int index) {
        if (atomList.size() == 0) {
            throw new IllegalArgumentException("empty list of Atoms");
        }
        int candidate = atomList.size()/2;
        int max = atomList.size();
        int min = 0;
        while (true) {
            Atom a = (Atom)atomList.get(candidate);
            int candidateIndex = a.getAddress();
            if (candidateIndex == index) {
                return a;
            }
            if (min == max - 1) {
                throw new IllegalArgumentException("could not find desired index in list of Atoms");
            }
            if (index > candidateIndex) {
                min = candidate;
                candidate = (max+min)/2;
                continue;
            }
            if (index < candidateIndex) {
                max = candidate;
                candidate = (max+min)/2;
                continue;
            }
            throw new IllegalArgumentException("you seem to have serialized multiple atoms with the same index but different signs");
        }
    }
}
