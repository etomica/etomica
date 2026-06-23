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

import java.util.Arrays;

/**
 * @author kofke
 *
 * Extension of Box that forms and holds a PairSet instance for all of the
 * atoms in the box.  Also instantiates box with a NONE boundary type.
 */
public class BoxCluster extends Box {

    protected boolean isTrial;
    protected CoordinatePairSet cPairSet, cPairTrialSet, cPairSetTmp;
    protected AtomPairSet aPairSet;
    protected long cPairID;
    protected final ClusterWeight sampleCluster;
    protected final Space space;
    protected IMoleculePositionDefinition positionDefinition;
    protected int[] types, idMap;
    protected MoleculeArrayList moleculesCluster;
    protected AtomArrayList atomsCluster;

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

    public IMoleculeList getMoleculeList() {
        // if we have an out-of-order mixture, return the molecules in the appropriate
        // order according to the diagram
        if (idMap == null) return super.getMoleculeList();
        if (moleculesCluster != null) return moleculesCluster;
        IMoleculeList simpleList = super.getMoleculeList();
        moleculesCluster = new MoleculeArrayList(simpleList.size());
        for (int i : idMap) {
            moleculesCluster.add(simpleList.get(i));
        }
        return moleculesCluster;
    }

    public IAtomList getLeafList() {
        // if we have an out-of-order mixture, return the atoms in the appropriate
        // order according to the diagram
        if (idMap == null) return super.getLeafList();
        if (atomsCluster != null) return atomsCluster;
        IAtomList simpleList = super.getLeafList();
        if (simpleList.size() != getMoleculeList().size()) {
            atomsCluster = new AtomArrayList(0);
            return atomsCluster;
        }
        atomsCluster = new AtomArrayList(simpleList.size());
        for (int i : idMap) {
            atomsCluster.add(simpleList.get(i));
        }
        return atomsCluster;
    }

    /**
     * Sets the type (species) index of each point in our diagram.
     * The alternate root point is not included, but is assumed to be
     * the same type as point 0.
     */
    public void setTypes(int[] types) {
        int nSpecies = moleculeLists.length;
        int nPoints = getMoleculeList().size()-1;
        // maps point index in the diagram to molecule index in box
        idMap = new int[nPoints+1];
        // maps given molecule from a given species to global index
        int[][] idMap1 = new int[nSpecies][nPoints+1];
        int k = 0;
        for (int i=0; i<nSpecies; i++) {
            int nti = moleculeLists[i].size();
            for (int j=0; j<nti; j++) {
                idMap1[i][j] = k;
                k++;
            }
        }
        int[] speciesUsed = new int[nSpecies];
        for (int i=0; i<nPoints; i++) {
            int ti = types[i];
            idMap[i] = idMap1[ti][speciesUsed[ti]];
            speciesUsed[ti]++;
        }
        // last diagram point is our alternate root
        idMap[nPoints] = idMap1[types[0]][speciesUsed[types[0]]];
        System.out.println("idMap: "+ Arrays.toString(idMap));
        // ensure we reconstruct our lists
        cPairSet = null;
        trialNotify();
        acceptNotify();
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
}
