/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.*;
import etomica.atom.IAtom;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;

public class BoxEventManager implements IBoxEventManager, java.io.Serializable {

    private transient final LinkedList<IBoxListener> intervalListeners = new LinkedList<IBoxListener>();
    private transient final ArrayList<Boolean> serial = new ArrayList<Boolean>();
    private final Box box;

    public BoxEventManager(Box _box) {
        box = _box;
    }
    
    public synchronized void addListener(IBoxListener newListener) {

        if(newListener == null) throw new NullPointerException("Cannot add null as a listener to Box");
        if (intervalListeners.contains(newListener)) {
            throw new RuntimeException(newListener+" is already an interval action");
        }
        intervalListeners.add(0, newListener);
        serial.add(0, true);
    }

    public synchronized void removeListener(IBoxListener listener) {
        intervalListeners.remove(listener);
    }

    public synchronized void moleculeAdded(IMolecule molecule) {
        IBoxMoleculeEvent event = new BoxMoleculeEvent(box, molecule);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).boxMoleculeAdded(event);
        }
    }
    
    public synchronized void moleculeRemoved(IMolecule molecule) {
        IBoxMoleculeEvent event = new BoxMoleculeEvent(box, molecule);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).boxMoleculeRemoved(event);
        }
    }

    public synchronized void globalAtomLeafIndexChanged(int index) {
        IBoxIndexEvent event = new BoxIndexEvent(box, index);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).boxGlobalAtomLeafIndexChanged(event);
        }
    }
    
    public synchronized void atomLeafIndexChanged(IAtom atom, int index) {
        IBoxAtomIndexEvent event = new BoxAtomIndexEvent(box, atom, index);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).boxAtomLeafIndexChanged(event);
        }
    }
    
    public synchronized void numberMolecules(ISpecies species, int count) {
        IBoxMoleculeCountEvent event = new BoxMoleculeCountEvent(box, species, count);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).boxNumberMolecules(event);
        }
    }
    
    public synchronized void moleculeIndexChanged(IMolecule molecule, int index) {
        IBoxMoleculeIndexEvent event = new BoxMoleculeIndexEvent(box, molecule, index);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).boxMoleculeIndexChanged(event);
        }
    }

    
    private void writeObject(java.io.ObjectOutputStream out)
    throws IOException
    {

        out.defaultWriteObject();
        
        // write # of listeners that will be serialized
        out.writeInt(intervalListeners.size());

        for(int i = 0; i < intervalListeners.size(); i++) {

            //skip transient listeners
            if (serial.get(i) == true) {
                out.writeObject(intervalListeners.get(i));
            }
        }
    }

    private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
    {

        in.defaultReadObject();
        
        // read the listener count
        int count = in.readInt();

        for (int i=0; i<count; i++) {
            addListener((IBoxListener)in.readObject());
        }
    }


}
