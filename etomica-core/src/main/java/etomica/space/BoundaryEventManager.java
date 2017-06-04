/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;

import etomica.api.IBoundaryEvent;
import etomica.api.IBoundaryEventManager;
import etomica.api.IBoundaryListener;

public class BoundaryEventManager implements IBoundaryEventManager, java.io.Serializable {

    private transient final LinkedList<IBoundaryListener> intervalListeners = new LinkedList<IBoundaryListener>();
    private transient final ArrayList<Boolean> serial = new ArrayList<Boolean>();
    
    public synchronized void addListener(IBoundaryListener newListener) {
        if(newListener == null) throw new NullPointerException("Cannot add null as a listener to Box");
        if (intervalListeners.contains(newListener)) {
            throw new RuntimeException(newListener+" is already an interval action");
        }
        intervalListeners.add(newListener);
        serial.add(0, true);
    }


    public synchronized void removeListener(IBoundaryListener listener) {
        intervalListeners.remove(listener);
    }

    public void inflate(Boundary boundary) {
        IBoundaryEvent event = new BoundaryEvent(boundary);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).boundaryInflate(event);
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
            addListener((IBoundaryListener)in.readObject());
        }
    }
}
