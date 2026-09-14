package etomica.models.water;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.units.Electron;

/**
 * Rigid four-site TIP4P/2005 water conformation.
 *
 * Atom ordering:
 *     H1, H2, O, M
 *
 * TIP4P/2005 geometry:
 *     O-H distance = 0.9572 Angstrom
 *     H-O-H angle  = 104.52 degrees
 *     O-M distance = 0.1546 Angstrom
 *
 * Charges:
 *     H1 = +0.5564 e
 *     H2 = +0.5564 e
 *     O  =  0.0000 e
 *     M  = -1.1128 e
 */
public class ConformationWaterTIP4P2005 implements IConformation, java.io.Serializable {

    private static final long serialVersionUID = 1L;

    protected final Space space;

    public static final double bondLengthOH = 0.9572;

    public static final double angleHOH = 104.52 * Math.PI / 180.0;

    public static final double rOM = 0.1546;

    public static final double[] Echarge = new double[4];

    static {
        ConformationWaterTIP4P2005.Echarge[SpeciesWater4P.indexH1] = Electron.UNIT.toSim(+0.5564);
        ConformationWaterTIP4P2005.Echarge[SpeciesWater4P.indexH2] = Electron.UNIT.toSim(+0.5564);
        ConformationWaterTIP4P2005.Echarge[SpeciesWater4P.indexO] = Electron.UNIT.toSim(0.0);
        ConformationWaterTIP4P2005.Echarge[SpeciesWater4P.indexM] = Electron.UNIT.toSim(-1.1128);
    }

    public ConformationWaterTIP4P2005(Space space) {
        this.space = space;
    }

    @Override
    public void initializePositions(IAtomList list) {

        if (list.size() != 4) {
            throw new IllegalArgumentException("TIP4P/2005 requires exactly four sites: " + "H1, H2, O, and M. Received " + list.size() + " sites.");
        }

        // Place oxygen at the origin.
        IAtom o = list.get(SpeciesWater4P.indexO);
        o.getPosition().E(new double[] {0.0, 0.0, 0.0});

        /*
         * Place the two hydrogen atoms symmetrically around
         * the positive y-axis. The y-axis is the H-O-H bisector.
         */
        double halfAngle = 0.5 * angleHOH;

        double x = bondLengthOH * Math.sin(halfAngle);

        double y = bondLengthOH * Math.cos(halfAngle);

        IAtom h1 = list.get(SpeciesWater4P.indexH1);
        h1.getPosition().E(new double[] {-x, y, 0.0});

        IAtom h2 = list.get(SpeciesWater4P.indexH2);
        h2.getPosition().E(new double[] {+x, y, 0.0});

        /*
         * The negative M site lies on the H-O-H bisector,
         * 0.1546 Angstrom from oxygen toward the hydrogens.
         */
        IAtom m = list.get(SpeciesWater4P.indexM);
        m.getPosition().E(new double[] {0.0, rOM, 0.0});
    }
}