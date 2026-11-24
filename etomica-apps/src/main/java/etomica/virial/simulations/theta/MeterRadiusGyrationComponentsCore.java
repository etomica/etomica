

package etomica.virial.simulations.theta;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
//import etomica.data.types.DataVector;
import etomica.data.types.DataDoubleArray;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Length;

/**
 * Core-aware version of MeterRadiusGyrationComponents. The stock version
 * weights every atom equally (weight=1) in both the center-of-mass and
 * the Rg tensor sum, exactly like the stock MeterRadiusGyration (confirmed
 * from source). This version instead follows Andrew's "replace the sum
 * over core atoms with the integral" approach: the core atom is weighted
 * by coreWeight (its effective atom count, M_core = rho*(4/3)*pi*s^3, per
 * the same convention already used for the scalar MeterRadiusGyrationCore),
 * and its own internal spread about its own center is added directly to
 * the tensor's diagonal.
 *
 * For a solid sphere, the scalar spread Rg^2=(3/5)*s^2 is isotropic, so it
 * splits evenly across the three diagonal components: each of xx, yy, zz
 * gets (1/5)*s^2 added (three components * (1/5)s^2 = (3/5)s^2, matching
 * the scalar result exactly). No off-diagonal (xy, xz, yz) contribution,
 * since a sphere has no preferred direction.
 *
 * coreWeight is NOT computed internally -- same reasoning as
 * MeterRadiusGyrationCore: whether "one effective core atom" should count
 * the same as one monomer bead in this weighted sum is a physical/modeling
 * decision, not something this class should assume silently.
 */
public class MeterRadiusGyrationComponentsCore implements IDataSource {
    private final Box box;
    private final AtomType coreType;
    private final double coreRadius;
    private final double coreWeight;
    private final Vector cm, realPos;
    private final Vector dr;
    //private final DataVector data;
    private final DataDoubleArray data;
    private final Tensor t;
    //private final DataVector.DataInfoVector datainfo;
    private final DataDoubleArray.DataInfoDoubleArray datainfo;
    private final DataTag tag;

    public MeterRadiusGyrationComponentsCore(Box box, AtomType coreType, double coreRadius, double coreWeight) {
        this.box = box;
        this.coreType = coreType;
        this.coreRadius = coreRadius;
        this.coreWeight = coreWeight;
        cm = box.getSpace().makeVector();
        realPos = box.getSpace().makeVector();
        dr = box.getSpace().makeVector();
        //data = new DataVector(box.getSpace());
        data = new DataDoubleArray(new int[]{4});
        t = box.getSpace().makeTensor();
        //datainfo = new DataVector.DataInfoVector("Radius Of Gyration (core-aware)", Length.DIMENSION, box.getSpace());
        datainfo = new DataDoubleArray.DataInfoDoubleArray("Radius of Gyration (core-aware)", Length.DIMENSION, new int[]{4});
        tag = new DataTag();
        datainfo.addTag(tag);
    }

    @Override
    public IData getData() {
        if (box == null) {
            throw new IllegalStateException("must call setBox before using meter");
        }
        t.E(0);
        Boundary boundary = box.getBoundary();
        double weightTot = 0.0;

        // (1/5)*s^2 -- the per-diagonal-component share of the core's own
        // internal Rg^2 = (3/5)*s^2 spread, weighted by coreWeight
        double coreDiagTerm = coreWeight * coreRadius * coreRadius / 5.0;

        for (IMolecule molecule : box.getMoleculeList()) {
            IAtomList childList = molecule.getChildList();
            if (childList.size() < 2) {
                continue;
            }

            // ---- pass 1: weighted center of mass ----
            IAtom firstAtom = childList.get(0);
            realPos.E(firstAtom.getPosition());
            double w0 = (firstAtom.getType() == coreType) ? coreWeight : 1.0;
            cm.E(realPos);
            cm.TE(w0);
            double wSum = w0;
            Vector prevPosition = firstAtom.getPosition();

            for (int iChild = 1; iChild < childList.size(); iChild++) {
                IAtom a = childList.get(iChild);
                Vector position = a.getPosition();
                dr.Ev1Mv2(position, prevPosition);
                boundary.nearestImage(dr);
                realPos.PE(dr);
                double w = (a.getType() == coreType) ? coreWeight : 1.0;
                Vector weightedPos = box.getSpace().makeVector();
                weightedPos.E(realPos);
                weightedPos.TE(w);
                cm.PE(weightedPos);
                wSum += w;
                prevPosition = position;
            }
            cm.TE(1.0 / wSum);

            // ---- pass 2: weighted tensor sum, plus core's own internal spread ----
            prevPosition = firstAtom.getPosition();
            realPos.E(firstAtom.getPosition());
            for (int iChild = 0; iChild < childList.size(); iChild++) {
                IAtom a = childList.get(iChild);
                Vector position = a.getPosition();
                dr.Ev1Mv2(position, prevPosition);
                boundary.nearestImage(dr);
                realPos.PE(dr);
                dr.Ev1Mv2(realPos, cm);

                boolean isCore = (a.getType() == coreType);
                double w = isCore ? coreWeight : 1.0;

                // t += w * (dr outer dr) -- built from PEv1v2 by pre-scaling dr,
                // since Tensor has no direct "scaled outer product add" method
                Vector scaledDr = box.getSpace().makeVector();
                scaledDr.E(dr);
                scaledDr.TE(w);
                t.PEv1v2(scaledDr, dr);

                if (isCore) {
                    // add (1/5)*s^2 (weighted) to each diagonal component only --
                    // isotropic spread, no off-diagonal contribution
                    for (int k = 0; k < t.D(); k++) {
                        t.PE(k, k, coreDiagTerm);
                    }
                }
                prevPosition = position;
            }
            weightTot += wSum;
        }
        t.TE(1.0 / weightTot);
        Matrix m = new Matrix(t.toArray(), t.D());
        EigenvalueDecomposition eigen = new EigenvalueDecomposition(m);
       // data.x.E(eigen.getRealEigenvalues());
        double[] ev = eigen.getRealEigenvalues();
        double[] out = data.getData();
        out[0] = ev[0];
        out[1] = ev[1];
        out[2] = ev[2];
        out[3] = ev[0] + ev[1] + ev[2];

        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return datainfo;
    }
}