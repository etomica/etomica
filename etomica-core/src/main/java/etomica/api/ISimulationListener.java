/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

public interface ISimulationListener {

    public void simulationBoxAdded(ISimulationBoxEvent e);
    
    public void simulationBoxRemoved(ISimulationBoxEvent e);
    
    public void simulationSpeciesAdded(ISimulationSpeciesEvent e);
    
    public void simulationSpeciesRemoved(ISimulationSpeciesEvent e);
    
    public void simulationSpeciesIndexChanged(ISimulationSpeciesIndexEvent e);
    
    public void simulationSpeciesMaxIndexChanged(ISimulationIndexEvent e);
    
    public void simulationAtomTypeIndexChanged(ISimulationAtomTypeIndexEvent e);

    public void simulationAtomTypeMaxIndexChanged(ISimulationIndexEvent e);
}
