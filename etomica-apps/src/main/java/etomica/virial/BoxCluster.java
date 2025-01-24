/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculeArrayList;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.virial.cluster.ClusterWeight;

import java.util.Map;

/**
 * @author kofke
 *
 * Extension of Box that forms and holds a PairSet instance for all of the
 * atoms in the box.  Also instantiates box with a NONE boundary type.
 */
public class BoxCluster extends Box {

	/**
	 * Constructor for BoxCluster.
	 */
	public BoxCluster(ClusterWeight cluster, Space _space) {
        this(cluster, _space, 0);
    }

    public BoxCluster(ClusterWeight cluster, Space _space, double L) {
        super(L == 0 ? new BoundaryRectangularNonperiodic(_space) : new BoundaryRectangularPeriodic(_space, L), _space);
        sampleCluster = cluster;
        this.space = _space;
	}

    /**
     * Set the mapping from molecule order in the box to molecule order in the diagram.
     * If you have species A and species B with diagram AABBA then the map should be
     *   0 => 0
     *   1 => 1
     *   2 => 3
     *   3 => 4
     *   4 => 2
     *
     *   This shouldn't be necessary for biconnected diagrams because we include
     *   all permutations, but can be useful for handling flexible correction
     *   diagrams.
     */
    public void setMoleculeMap(Map<Integer,Integer> idMap) {
        this.idMap = idMap;
    }

	public void setPositionDefinition(IMoleculePositionDefinition positionDefinition){
	    this.positionDefinition = positionDefinition;
        if(cPairSet != null && cPairSet instanceof CoordinatePairMoleculeSet ){
            ((CoordinatePairMoleculeSet) cPairSet).setPositionDefinition(positionDefinition);
            ((CoordinatePairMoleculeSet) cPairTrialSet).setPositionDefinition(positionDefinition);
        }
    }

    /**
     * returns the current coordinate pair set
     */
	public CoordinatePairSet getCPairSet() {
		return isTrial ? cPairTrialSet : cPairSet;
	}

	public long getCPairID() {
	    return cPairID;
	}

    public AtomPairSet getAPairSet() {
        return aPairSet;
    }

	/**
     * returns the cluster used for sampling in this box
	 */
    public ClusterWeight getSampleCluster() {
        return sampleCluster;
    }

    /**
     * Inform the box that a trial move has been made so it can update
     * the coordinate pairs.  If molecule is not null, only coordinate pairs 
     * containing that atom are updated.
     */
    public void trialNotify() {
        // atom(s) have been moved.  leave cPairSet as is and update
        // cPairTrialSet and set a flag to use it.
        isTrial = true;
        // increase ID to notify clusters to recalculate value
        if(cPairSet == null) {
            // assume 1 species
            IMoleculeList molecules = getMoleculeList();
            IAtomList leafAtoms = getLeafList();
            boolean doAtoms = molecules.size() == leafAtoms.size();
            if (idMap != null) {
                MoleculeArrayList mNew = new MoleculeArrayList();
                AtomArrayList aNew = new AtomArrayList();
                for (int i=0; i<molecules.size(); i++) {
                    int j = idMap.get(i);
                    mNew.add(molecules.get(j));
                    if (doAtoms) {
                        aNew.add(leafAtoms.get(j));
                    }
                }
                molecules = mNew;
                leafAtoms = aNew;
            }
            if (doAtoms) {
                cPairSet = new CoordinatePairLeafSet(leafAtoms,space);
                cPairTrialSet = new CoordinatePairLeafSet(leafAtoms,space);
            }
            else {
                cPairSet = new CoordinatePairMoleculeSet(molecules,space);
                if(positionDefinition!=null)((CoordinatePairMoleculeSet) cPairSet).setPositionDefinition(positionDefinition);
                cPairTrialSet = new CoordinatePairMoleculeSet(molecules,space);
                if(positionDefinition!=null)((CoordinatePairMoleculeSet) cPairTrialSet).setPositionDefinition(positionDefinition);
            }
            aPairSet = new AtomPairSet(molecules);
        }

        cPairID++;
        cPairTrialSet.reset(cPairID);
    }

    /**
     * Informs the box that the trial was accepted so it will keep the new 
     * coordinate pairs.
     */
	public void acceptNotify() {
        if (!isTrial) {
            throw new IllegalStateException("you weren't in a trial!");
        }
	    // move was accepted.  swap out trial cPairSet and cPairTrialSet since
		// cPairTrialSet is already up-to-date
		isTrial = false;
		cPairSetTmp = cPairSet;
		cPairSet = cPairTrialSet;
		cPairTrialSet = cPairSetTmp;
	}

    /**
     * Informs the box that the trial was rejected so it will go back to 
     * the old coordinate pairs.
     */
	public void rejectNotify() {
        if (!isTrial) {
            throw new IllegalStateException("you weren't in a trial!");
        }
		// move was rejected.  stop using cPairTrialSet.
		isTrial = false;
	}

 	protected boolean isTrial;
	protected CoordinatePairSet cPairSet, cPairTrialSet, cPairSetTmp;
    protected AtomPairSet aPairSet;
    protected long cPairID;
	protected final ClusterWeight sampleCluster;
	protected final Space space;
	protected IMoleculePositionDefinition positionDefinition;
    protected Map<Integer,Integer> idMap;
}
