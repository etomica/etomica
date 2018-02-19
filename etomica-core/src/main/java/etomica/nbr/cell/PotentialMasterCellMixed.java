/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell;

import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculation;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * A PotentialMaster for use with mixtures having very different iteraction ranges.
 * Short interactions can be handled by cells while longer-range interactions are
 * handled with brute force -- all pairs of atoms are computed.  The cells can then
 * be sized to optimize the short interactions.
 *
 * @author Andrew Schultz
 */
public class PotentialMasterCellMixed extends PotentialMasterCell {

    protected final Set<ISpecies> cellSpecies;
    protected Box box;
    protected Map<AtomType, Map<AtomType, IPotentialAtomic>> unrangedPotentials;
    protected final AtomPair pair;
    protected IteratorDirective.Direction direction;

    public PotentialMasterCellMixed(Simulation sim, double range, Space _space) {
        this(sim, range, new BoxAgentSourceCellManager(sim, null, _space), _space);
    }

    public PotentialMasterCellMixed(Simulation sim, double range,
                                    BoxAgentSourceCellManager boxAgentSource, Space _space) {
        this(sim, range, boxAgentSource, new BoxAgentManager<NeighborCellManager>(boxAgentSource, NeighborCellManager.class), _space);
    }

    public PotentialMasterCellMixed(Simulation sim, double range, BoxAgentSourceCellManager boxAgentSource,
                                    BoxAgentManager<NeighborCellManager> agentManager, Space _space) {
        super(sim, range, boxAgentSource, agentManager, _space);
        cellSpecies = new HashSet<>();
        unrangedPotentials = new HashMap<>();
        pair = new AtomPair();
        setRange(range);
    }

    public void setHandledByCells(ISpecies species, boolean isHandled) {
        if (isHandled) cellSpecies.add(species);
        else if (cellSpecies.contains(species)) cellSpecies.remove(species);
    }

    public void addUnrangedPotential(AtomType type1, AtomType type2, IPotentialAtomic p2) {
        Map<AtomType, IPotentialAtomic> pMap = unrangedPotentials.get(type1);
        if (pMap == null) {
            pMap = new HashMap<>();
            unrangedPotentials.put(type1, pMap);
        }
        pMap.put(type2, p2);

        pMap = unrangedPotentials.get(type2);
        if (pMap == null) {
            pMap = new HashMap<>();
            unrangedPotentials.put(type2, pMap);
        }
        pMap.put(type1, p2);

    }

    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        this.box = box;
        this.direction = id.direction();
        for (Map<AtomType, IPotentialAtomic> pMap : unrangedPotentials.values()) {
            for (IPotentialAtomic p : pMap.values()) {
                p.setBox(box);
            }
        }
        super.calculate(box, id, pc);
    }

    protected void calculate(IAtom atom, PotentialCalculation pc) {
        boolean atomCells = cellSpecies.contains(atom.getType().getSpecies());
        Map<AtomType, IPotentialAtomic> pMap = unrangedPotentials.get(atom.getType());
        if (atomCells) {
            // some of the atoms interactions are handled with cell lists
            // do that now
            super.calculate(atom, pc);
            if (pMap == null) return;
            pair.atom0 = atom;
            for (int i = 0; i < simulation.getSpeciesCount(); i++) {
                ISpecies s = simulation.getSpecies(i);
                if (cellSpecies.contains(s)) {
                    // we already handled these interactions with the cell lists
                    continue;
                }
                IMoleculeList molecules = box.getMoleculeList(s);
                for (int j = 0; j < molecules.getMoleculeCount(); j++) {
                    IAtomList atoms = molecules.getMolecule(j).getChildList();
                    for (int k = 0; k < atoms.getAtomCount(); k++) {
                        IAtom kAtom = atoms.getAtom(k);
                        IPotentialAtomic p2 = pMap.get(kAtom.getType());
                        if (p2 == null) continue;
                        pair.atom1 = kAtom;
//                        System.out.println("brute-force handling "+pair);
                        pc.doCalculation(pair, p2);
                    }
                }
            }
        } else {
            if (pMap == null) return;
            // not handled by cells at all.  loop over all atoms
            IAtomList atoms = box.getLeafList();
            pair.atom0 = atom;
            int i = atom.getLeafIndex();
            int start = direction != IteratorDirective.Direction.UP ? 0 : i + 1;
            int stop = direction != IteratorDirective.Direction.DOWN ? atoms.getAtomCount() : i;
            for (int j = start; j < stop; j++) {
                if (j == i) continue;
                IAtom jAtom = atoms.getAtom(j);
                IPotentialAtomic p2 = pMap.get(jAtom.getType());
                if (p2 == null) continue;
                pair.atom1 = jAtom;
//                System.out.println("brute-force handling "+pair);
                pc.doCalculation(pair, p2);
            }
        }
    }

}
