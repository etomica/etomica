package etomica.virial;

public class ClusterWheatleySoftDerivativesMix extends ClusterWheatleySoftDerivatives {

    protected final MayerFunction[][] mixF;
    protected final int[] nTypes;
    protected final MayerFunction[][] fMap;

    public ClusterWheatleySoftDerivativesMix(int nPoints, int[] nTypes, MayerFunction[][] f, double tol, int nDer){
        super(nPoints,null,tol,nDer);
        this.nTypes=nTypes;
        mixF = f;
        fMap = new MayerFunction[nPoints][nPoints];
        if(tol!=0){setTolerance(tol);}
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

    public void setTolerance(double newTol) {
        if(nTypes==null)return;
        if(newTol!=0){
            clusterBD = new ClusterWheatleySoftDerivativesMixBD(n, nTypes, mixF, -3*(int)Math.log10(newTol),nDer);
            clusterBD.setDoCaching(false);
            clusterBD.setPrecisionLimit(300);
        }
        else{
            clusterBD = null;
        }
        tol = newTol;
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleySoftDerivativesMix c = new ClusterWheatleySoftDerivativesMix(n, nTypes, mixF, tol,nDer);
        c.setTemperature(1/beta);
        c.setDoCaching(doCaching);
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
                if (debug && (Double.isNaN(ff) || Double.isInfinite(ff))) {
                    System.err.println("oops in updateF "+i+" "+j+" "+ff);
                }
//                if (Math.abs(ff) < 1e-14) ff = 0;
                fQ[(1<<i)|(1<<j)][0] = ff+1;
            }
        }
    }
}
