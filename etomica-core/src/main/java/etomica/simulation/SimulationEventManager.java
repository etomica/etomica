/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import java.io.IOException;
import java.util.LinkedList;

import etomica.atom.IAtomType;
import etomica.box.Box;
import etomica.api.ISimulationAtomTypeIndexEvent;
import etomica.api.ISimulationBoxEvent;
import etomica.api.ISimulationEventManager;
import etomica.api.ISimulationIndexEvent;
import etomica.api.ISimulationListener;
import etomica.api.ISimulationSpeciesEvent;
import etomica.api.ISimulationSpeciesIndexEvent;
import etomica.api.ISpecies;

public class SimulationEventManager implements ISimulationEventManager {

    private transient final LinkedList<ISimulationListener> intervalListeners = new LinkedList<ISimulationListener>();
    private final Simulation simulation;
    
    public SimulationEventManager(Simulation sim) {
        simulation = sim;
    }

    /* (non-Javadoc)
     * @see etomica.util.IEventManager#addListener(java.lang.Object)
     */
    public synchronized void addListener(ISimulationListener listener) {
        if(listener == null) throw new NullPointerException("Cannot add null as a listener to Box");
//        if (intervalListeners.contains(listener)) {
//            throw new RuntimeException(listener+" is already an interval action");
//        }
        intervalListeners.add(listener);
    }

    public synchronized void boxAdded(Box box) {
        ISimulationBoxEvent e = new SimulationBoxEvent(simulation, box);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationBoxAdded(e);
        }
    }
    
    public synchronized void boxRemoved(Box box) {
        ISimulationBoxEvent e = new SimulationBoxEvent(simulation, box);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationBoxRemoved(e);
        }
    }
    
    public synchronized void speciesAdded(ISpecies species) {
        ISimulationSpeciesEvent e = new SimulationSpeciesEvent(simulation, species);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationSpeciesAdded(e);
        }
    }
    
    public synchronized void speciesRemoved(ISpecies species) {
        ISimulationSpeciesEvent e = new SimulationSpeciesEvent(simulation, species);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationSpeciesRemoved(e);
        }
    }
    
    public synchronized void speciesIndexChanged(ISpecies species, int index) {
        ISimulationSpeciesIndexEvent e = new SimulationSpeciesIndexEvent(simulation, species, index);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationSpeciesIndexChanged(e);
        }
    }
    
    public synchronized void speciesMaxIndexChanged(int index) {
        ISimulationIndexEvent e = new SimulationIndexEvent(simulation, index);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationSpeciesMaxIndexChanged(e);
        }
    }
    
    public synchronized void atomTypeIndexChanged(IAtomType atomType, int index) {
        ISimulationAtomTypeIndexEvent e = new SimulationAtomTypeIndexEvent(simulation, atomType, index);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationAtomTypeIndexChanged(e);
        }
    }

    public synchronized void atomTypeMaxIndexChanged(int index) {
        ISimulationIndexEvent e = new SimulationIndexEvent(simulation, index);
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationAtomTypeMaxIndexChanged(e);
        }
    }
    
    /* (non-Javadoc)
     * @see etomica.util.IEventManager#removeListener(java.lang.Object)
     */
    public synchronized void removeListener(ISimulationListener listener) {
        intervalListeners.remove(listener);
    }
    
    private void writeObject(java.io.ObjectOutputStream out)
    throws IOException
    {
        
        out.defaultWriteObject();
        
        // write # of listeners that will be serialized
        out.writeInt(intervalListeners.size());

        for(int i = 0; i < intervalListeners.size(); i++) {

            out.writeObject(intervalListeners.get(i));

        }
        
        
    }

    private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
    {
        in.defaultReadObject();
        
        // read the listener count
        int count = in.readInt();

        for (int i=0; i<count; i++) {
            addListener((ISimulationListener)in.readObject());
        }
    }

}
