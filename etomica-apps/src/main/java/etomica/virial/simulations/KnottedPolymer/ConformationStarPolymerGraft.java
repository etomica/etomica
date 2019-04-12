package etomica.virial.simulations.KnottedPolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * General Class for the conformation of Star polymers
 */
public class ConformationStarPolymerGraft implements IConformation {
    protected final Space space;
    public int f, l;
    protected double sigma = 1.0;
    protected Vector[] graftVectors;

    public ConformationStarPolymerGraft(Space space, int f, int l) {
        this.space = space;
        this.f = f;
        this.l = l;
    }

    @Override
    public void initializePositions(IAtomList atomList) {
        IAtom core = atomList.getAtoms().get(0);
        core.getPosition().E(0);
        Vector corePos = core.getPosition();
        graftVectors = new Vector[12];

        if (f <= 12) {
            for (int i = 0; i < 12; i++) {
                graftVectors[i] = space.makeVector();
            }
            double phi = (1 + Math.sqrt(4)) / 2;
            double x1 = 1 / Math.sqrt(phi * phi + 1);
            double x2 = phi * x1;
            graftVectors[0].E(new double[]{0, x1, x2});
            graftVectors[1].E(new double[]{0, x1, -x2});
            graftVectors[2].E(new double[]{0, -x1, x2});
            graftVectors[3].E(new double[]{0, -x1, -x2});
            graftVectors[4].E(new double[]{x1, x2, 0});
            graftVectors[5].E(new double[]{x1, -x2, 0});
            graftVectors[6].E(new double[]{-x1, x2, 0});
            graftVectors[7].E(new double[]{-x1, -x2, 0});
            graftVectors[8].E(new double[]{x2, 0, x1});
            graftVectors[9].E(new double[]{x2, 0, -x1});
            graftVectors[10].E(new double[]{-x2, 0, x1});
            graftVectors[11].E(new double[]{-x2, 0, -x1});
        }

        int iBead = 1;
        Vector dr = space.makeVector();
        Vector lastPos = space.makeVector();

        for (int j = 0; j < f; j++) {
            lastPos.E(corePos);
//            lastPos.PEa1Tv1(sigmaBead + bond, graftVectors[j]);
            for (int k = 0; k < l; k++) {
                IAtom atom = atomList.get(iBead);
                dr.E(graftVectors[j]);

                Vector p = atom.getPosition();
                p.E(lastPos);
                p.PEa1Tv1((k == 0 ? 1.15 : 1.15) * sigma, dr);
                iBead++;
                lastPos.E(p);
            }
        }

    }
}
