package etomica.modules.interfacial;

import java.util.Random;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.potential.PotentialSoft;
import etomica.space.Tensor;

public class P1Lrc implements PotentialSoft {

    public P1Lrc(PotentialCalculationForcePressureBinSum forceSum) {
        this.forceSum = forceSum;
        gradient = new IVector[1];
        gradient[0] = forceSum.space.makeVector();
        random = new Random();
    }
    
    public IVector[] gradient(IAtomList atoms) {
        return gradient(atoms, null);
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        double binSize = forceSum.binSize;
        double[] phi = forceSum.fLrc;
        double[] xpi = forceSum.xVirialLrc;
        double[] yzpi = forceSum.yzVirialLrc;
        double[] densityProfile = forceSum.densityProfile;
        IVector dim = box.getBoundary().getDimensions();
        double L = dim.x(0);
        IAtomPositioned a = (IAtomPositioned)atoms.getAtom(0);
        double x = a.getPosition().x(0);
        x -= Math.round(x/L) * L;
        int iBin = (int) ((x + 0.5 * L) / binSize);
        int nBins = densityProfile.length;
        int halfNBins = nBins/2;
        double fx = 0;
        double xVirial = 0;
        double yzVirial = 0;
        for (int i=0; i<nBins; i++) {
            int di = iBin - i;
            int mdi = di;
            if (di < 0) mdi = -mdi;
            if (mdi == halfNBins) continue;
            if (mdi > halfNBins) {
                mdi = nBins - mdi;
                if (di > 0) {
                    di = -mdi;
                }
                else {
                    di = mdi;
                }
            }
            fx += di * binSize * binSize * densityProfile[i] * phi[mdi];
            xVirial += binSize * densityProfile[i] * xpi[mdi];
            yzVirial += binSize * densityProfile[i] * yzpi[mdi];
        }
        if (Double.isNaN(fx) || Double.isInfinite(fx)) {
            throw new RuntimeException("oops "+fx);
        }
        gradient[0].setX(0, -fx);
        if (pressureTensor != null) {
            pressureTensor.setComponent(0, 0, xVirial);
            pressureTensor.setComponent(1, 1, yzVirial);
            if (pressureTensor.D() == 3) {
                pressureTensor.setComponent(2, 2, yzVirial);
            }
        }
        if (random.nextDouble() > 1) {
            System.out.println(x+" "+fx+" "+xVirial+" "+yzVirial);
        }
        return gradient;
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    public double energy(IAtomList atoms) {
        double binSize = forceSum.binSize;
        double[] w = forceSum.uLrc;
        double[] densityProfile = forceSum.densityProfile;
        IVector dim = box.getBoundary().getDimensions();
        double L = dim.x(0);
        int nBins = densityProfile.length;
        if (Math.abs(L - binSize * nBins) > 1) return 0;
        IAtomPositioned a = (IAtomPositioned)atoms.getAtom(0);
        double x = a.getPosition().x(0);
        x -= Math.round(x/L) * L;
        int iBin = (int) ((x + 0.5 * L) / binSize);
        int halfNBins = nBins/2;
        double u = 0;
        for (int i=0; i<nBins; i++) {
            int di = i - iBin;
            if (di < 0) di = -di;
            if (di == halfNBins) continue;
            if (di > halfNBins) {
                di = nBins - di;
            }
            u += densityProfile[i] * w[di];
        }
        return u * binSize;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public int nBody() {
        return 1;
    }

    public void setBox(IBox box) {
        this.box = box;
    }

    protected final PotentialCalculationForcePressureBinSum forceSum;
    protected IBox box;
    protected final IVector[] gradient;
    protected final Random random;
}
