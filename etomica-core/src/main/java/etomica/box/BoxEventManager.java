/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.species.ISpecies;

import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

public class BoxEventManager {

    private final List<BoxEventListener> listeners = new CopyOnWriteArrayList<>();
    private final Box box;

    public BoxEventManager(Box box) {
        this.box = box;
    }

    public void addListener(BoxEventListener newListener) {

        if (newListener == null) throw new NullPointerException("Cannot add null as a listener to Box");
        if (listeners.contains(newListener)) {
            throw new RuntimeException(newListener + " is already an interval action");
        }
        listeners.add(newListener);
    }

    public void removeListener(BoxEventListener listener) {
        listeners.remove(listener);
    }

    public void moleculeAdded(IMolecule molecule) {
        BoxMoleculeEvent event = new BoxMoleculeEvent(box, molecule);
        for (BoxEventListener listener : listeners) {
            listener.boxMoleculeAdded(event);
        }
    }

    public void moleculeRemoved(IMolecule molecule) {
        BoxMoleculeEvent event = new BoxMoleculeEvent(box, molecule);
        for (BoxEventListener listener : listeners) {
            listener.boxMoleculeRemoved(event);
        }
    }

    public void globalAtomLeafIndexChanged(int index) {
        BoxIndexEvent event = new BoxIndexEvent(box, index);
        for (BoxEventListener listener : listeners) {
            listener.boxGlobalAtomLeafIndexChanged(event);
        }
    }

    public void atomLeafIndexChanged(IAtom atom, int index) {
        BoxAtomIndexEvent event = new BoxAtomIndexEvent(box, atom, index);
        for (BoxEventListener listener : listeners) {
            listener.boxAtomLeafIndexChanged(event);
        }
    }

    public void numberMolecules(ISpecies species, int count) {
        BoxMoleculeCountEvent event = new BoxMoleculeCountEvent(box, species, count);
        for (BoxEventListener listener : listeners) {
            listener.boxNumberMolecules(event);
        }
    }

    public void moleculeIndexChanged(IMolecule molecule, int index) {
        BoxMoleculeIndexEvent event = new BoxMoleculeIndexEvent(box, molecule, index);
        for (BoxEventListener listener : listeners) {
            listener.boxMoleculeIndexChanged(event);
        }
    }
}
