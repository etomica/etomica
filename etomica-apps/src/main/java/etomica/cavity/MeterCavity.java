package etomica.cavity;

import etomica.data.meter.MeterRDF;
import etomica.space.Boundary;
import etomica.space.Space;

/**
 * Meter collects RDF but is hardcoded to look at only those atoms that are
 * overlapping.
 */
public class MeterCavity extends MeterRDF {

    protected final P2HardSphereCavity p2;

    public MeterCavity(Space space, P2HardSphereCavity p2) {
        super(space, true);
        this.p2 = p2;
    }

    public void actionPerformed() {
        if (p2.pairedAtom1 == null) return;
        if (rData != xDataSource.getData() ||
                data.getLength() != rData.getLength() ||
                xDataSource.getXMax() != xMax) {
            reset();
        }

        double xMaxSquared = xMax * xMax;
        iterator.setBox(box);
        iterator.reset();
        Boundary boundary = box.getBoundary();
        // iterate over all pairs
        dr.Ev1Mv2(p2.pairedAtom1.getPosition(), p2.pairedAtom2.getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 < xMaxSquared) {
            int index = xDataSource.getIndex(Math.sqrt(r2));
            gSum[index]++;
        }
        callCount++;
    }
}
