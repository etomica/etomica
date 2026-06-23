package etomica.models.traPPE;

import etomica.molecule.IMolecule;
import etomica.potential.SiteReconstructor;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.virial.mcmove.MCMoveClusterAngleGeneral;

public class MCMoveClusterAngleGeneralImplicit extends MCMoveClusterAngleGeneral {

    protected final SiteReconstructor reconstructor;

    public MCMoveClusterAngleGeneralImplicit(PotentialCompute potentialCompute, Space space, ISpecies species, IntArrayList[] bonding, boolean oneSide, int[][][] triplets, IRandom random, double stepSize, SiteReconstructor reconstructor) {
        super(potentialCompute, space, species, bonding, oneSide, triplets, random, stepSize);
        this.reconstructor = reconstructor;
    }

    @Override
    public boolean doTrial() {
        boolean trial = super.doTrial();
        if (!trial) return false;
        IMolecule molecule = box.getMoleculeList().get(iMolecule);
        reconstructor.reconstructSites(molecule);
        return true;
    }
}
