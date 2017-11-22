/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.site;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.BoxCellManager;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;
import etomica.space.Space;

/**
 * Neighbor manager for system in which there is an unchanging, one-to-one
 * correspondence between atoms and lattice cells. Each leaf atom is assigned to
 * its own cell at the beginning of the simulation, and the cell-atom assignment
 * doesn't change through the course of the simulation. Atom neighbors are
 * determined by the structure and neighbor definition associated with the cell
 * lattice, and do not change through the simulation.
 * <p>
 * The expected use for this neighbor structure is in modeling of
 * lattice systems such as the Ising model in which the coordinate represents
 * something other than the spatial position of the atom; can also be applied to
 * diffusionless crystalline solids, such as the valence-force field model,
 * single-occupancy crystals, and perhaps high-density unconstrained crystals.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on May 23, 2005 by kofke
 */
public class NeighborSiteManager implements BoxCellManager, AgentSource<AtomSite> {

    /**
     * Constructs manager for neighbor cells in the given box. The number of
     * cells in each dimension is given by nCells. Position definition for each
     * atom is that given by its type (it is set to null in this class).
     */
    public NeighborSiteManager(final Box box, int nCells, Space _space) {
        space = _space;

        lattice = new CellLattice(space, box.getBoundary().getBoxSize(),
                AtomSite.FACTORY);
        int[] size = new int[space.D()];
        for (int i = 0; i < space.D(); i++)
            size[i] = nCells;
        lattice.setSize(size);
        siteIterator = new RectangularLattice.Iterator(space.D());
        siteIterator.setLattice(lattice);
        siteIterator.reset();

        agentManager = new AtomLeafAgentManager<AtomSite>(this,box);
    }

    /**
     * Returns the lattice used to define the neighbor structure.
     */
    public CellLattice getLattice() {
        return lattice;
    }

    /**
     * Should not be called, because cell assignments are made as atoms are created.
     * 
     * @throws RuntimeException
     *             if invoked
     */
    public void assignCellAll() {
        throw new RuntimeException("Cell assignments are made automagically.  This method isn't here.  Really.");
    }

    public AtomSite getSite(IAtom atom) {
        return agentManager.getAgent(atom);
    }

    public AtomSite makeAgent(IAtom atom, Box agentBox) {
        AtomSite site = (AtomSite)siteIterator.next();
        site.setAtom(atom);
        return site;
    }
    
    public void releaseAgent(AtomSite agent, IAtom atom, Box agentBox) {}

    private final CellLattice lattice;
    private final Space space;
    private final RectangularLattice.Iterator siteIterator;
    private final AtomLeafAgentManager<AtomSite> agentManager;
}
