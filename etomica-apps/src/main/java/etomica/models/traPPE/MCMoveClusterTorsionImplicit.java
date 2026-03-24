package etomica.models.traPPE;

import etomica.molecule.IMolecule;
import etomica.potential.SiteReconstructor;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.virial.mcmove.MCMoveClusterTorsion;

public class MCMoveClusterTorsionImplicit extends MCMoveClusterTorsion {
    protected final SiteReconstructor reconstructor;

    public MCMoveClusterTorsionImplicit(PotentialCompute potentialCompute, Space space, ISpecies species, IntArrayList[] bonding, int[][][] quads, IRandom random, double stepSize, SiteReconstructor reconstructor) {
        super(potentialCompute, space, species, bonding, quads, random, stepSize);
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
