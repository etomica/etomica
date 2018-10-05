package etomica.virial;

import java.math.BigDecimal;

public class ClusterWheatleySoftDerivativesMixBD extends  ClusterWheatleySoftDerivativesBD{

    protected final MayerFunction[][] mixF;
    protected final int[] nTypes;
    protected final MayerFunction[][] fMap;

    public ClusterWheatleySoftDerivativesMixBD(int nPoints, int[] nTypes, MayerFunction[][] f, int precision, int nDer) {
        super(nPoints, null, precision, nDer);
        this.nTypes = nTypes;
        mixF = f;
        fMap = new MayerFunction[nPoints][nPoints];
        int iType = 0, jType = 0;
        int iSum = nTypes[0], jSum = 0;
        for (int i=0; i<nPoints; i++) {
            while (i>=iSum) {
                iType++;
                iSum += nTypes[iType];
            }
            jSum = iSum;
            jType = iType;
            for (int j=i+1; j<nPoints; j++) {
                while (j>=jSum) {
                    jType++;
                    jSum += nTypes[jType];
                }
                fMap[i][j] = f[iType][jType];
            }
        }

    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleySoftDerivativesMixBD c = new ClusterWheatleySoftDerivativesMixBD(n, nTypes, mixF, mc.getPrecision(),nDer);
        c.setTemperature(1/beta);
        return c;
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();
        for (int i=0; i<mixF.length; i++) {
            for (int j=0; j<mixF[i].length; j++) {
                mixF[i][j].setBox(box);
            }
        }

        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                double ff = fMap[i][j].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
//                if (Math.abs(ff) < 1e-14) ff = 0;
                fQ[(1<<i)|(1<<j)][0] = new BigDecimal(ff).add(BDONE, mc);
            }
        }
    }

}

