package etomica.virial;



/**
 * This class calculates the sum of all chain clusters. Can be configured at construction to do rings instead.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterChain implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;
    protected final boolean doRing;
    
    protected final double[][] nC;
    protected final double[][] fValues;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    
    public ClusterChain(int nPoints, MayerFunction f) {
        this(nPoints, f, false);
    }
    
    public ClusterChain(int nPoints, MayerFunction f, boolean doRing) {
        this.n = nPoints;
        this.f = f;
        this.doRing = doRing;
        nf = 1<<n;  // 2^n

        nC = new double[n][nf];
        fValues = new double[n][n];
    }

    public ClusterAbstract makeCopy() {
        ClusterChain c = new ClusterChain(n, f);
        c.setTemperature(1/beta);
        return c;
    }
    
    public int pointCount() {
        return n;
    }


    public double value(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        int thisCPairID = cPairs.getID();
        if (thisCPairID == cPairID) {
            return value;
        }
        if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            value = lastValue;
            return value;
        }

        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;
        
        updateF(box);
        
        calcValue();
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            updateF(box);
            calcValue();
            throw new RuntimeException("oops");
        }
        return value;
    }

    public long numDiagrams() {
        double savedValue = value;
        for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                 fValues[i][j] = 1;
                 fValues[j][i] = 1;
            }
        }
        
        calcValue();
        long num = (int)Math.round(value);
        value = savedValue;
        return num;
    }


    /*
     * Computation of sum of chain or ring diagrams.
     */
    protected void calcValue() {
        
        //nC[m][i] is the number of chains beginning at 0 and ending at m, traversing all points in i
        
        //Start with all pairwise paths from 0 to each vertex
        for(int m=1; m<n; m++) {
            nC[m][(1<<m)|1] = fValues[m][0];
        }
        
        //All other paths
        //(could probably reduce memory by not including redundant first bit)
        for(int i=3; i<nf-1; i+=2) {//1-bit is always nonzero in i
            for(int m=1; m<n; m++) {//loop over indices not in i
                int im = 1<<m;
                if((im & i) != 0) continue;//skip if m is in i
                int index = i|im;
                nC[m][index] = 0;
                for(int k=1; k<n; k++) {//loop over indices in i
                    int ik = 1<<k;
                    if(ik > i) break;
                    if((ik & i) == 0) continue;//skip if k is not in i
                    nC[m][index] += fValues[m][k]*nC[k][i];
                }
            }
        }
        
        value = 0;
        if(doRing) {
            for(int m=1; m<n; m++) {
                value += nC[m][nf-1] * fValues[m][0];
            }
            
        } else { //chains
        
            //Sum chains in which first vertex is not a leaf.
            //Consider all partitions, counting paths beginning in one partition and ending in its complement
            for(int iS=3; iS<nf; iS+=4) {//keep 1 and 2 in i-partition to prevent double counting
                int iSComp = (nf-1)^iS;
                for(int m=1; m<n; m++) {
                    if(((1<<m)&iS) == 0) continue;//skip if m is not in iS
                    for(int k=2; k<n; k++) {
                        if(((1<<k)&iSComp) == 0) continue;//skip if k is not in iSComp
                        value += nC[m][iS] * nC[k][iSComp|1];
                    }
                }
            }
            
            //Sum chains where first vertex is a leaf
            for(int m=1; m<n; m++) {
                value += nC[m][nf-1];
            }
        
        }
        
    }
        
    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                fValues[i][j] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                fValues[j][i] = fValues[i][j];
            }
        }
    }
        
    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
    
    /**
     * Normalized probability for the separation of the centers of two spheres located at the
     * end of a chain of (n+1) overlapping spheres. All spheres of unit diameter.
     * @param n number of bonds in the chain, such that number of spheres is n+1
     * @param r end-to-end separation
     * @return value of normalized probability density
     */
    public double separationProbability(int n, double r) {

        if(r >= n) return 0;

        switch(n) {
            case 2:
                return Power(-2 + r,2) * (4 + r)/12.;
            case 3:
                if(r < 1) return (-525 + 315*Power(r,2) - 63*Power(r,4) + Power(r,6))/(9.*(-92 + 27*Math.log(3.)));
                else return -(Power(-3 + r,4)*(-6 + 27*r + 12*Power(r,2) + Power(r,3)))/(18.*r*(-92 + 27*Math.log(3.)));
            case 4:
                if(r < 2) return (-348160 + 184320*Power(r,2) - 56448*Power(r,4) + 15120*Power(r,5) + 960*Power(r,6) - 
                        540*Power(r,7) + 3*Power(r,9))/(4608.*(-179 + 256*ArcCoth(3)));
                else return -(Power(-4 + r,6)*(-144 + 224*r + 156*Power(r,2) + 24*Power(r,3) + Power(r,4)))/
                (4608.*r*(-179 + 256*ArcCoth(3)));
            case 5:
                if(r < 1) return -(146392675 - 63050130*Power(r,2) + 12657645*Power(r,4) - 1501500*Power(r,6) + 
                        96525*Power(r,8) - 1170*Power(r,10) + 3*Power(r,12))/
                        (180.*(-2773712 + 6640625*ArcCoth(4) + 3727*Math.log(3)));
                else if(r < 3) return (335430 - 75822175*r + 8925150*Power(r,2) + 14401530*Power(r,3) + 
                        20091500*Power(r,4) - 20765745*Power(r,5) + 5791500*Power(r,6) - 
                        42900*Power(r,7) - 244530*Power(r,8) + 32175*Power(r,9) + 1430*Power(r,10) - 
                        390*Power(r,11) + Power(r,13))/
                      (90.*r*(-2773712 + 6640625*ArcCoth(4) + 3727*Math.log(3)));
                else return -(Power(-5 + r,8)*(-3060 + 1825*r + 2260*Power(r,2) + 510*Power(r,3) + 40*Power(r,4) + 
                        Power(r,5)))/(360.*r*(-2773712 + 6640625*ArcCoth(4) + 3727*Math.log(3)));
            case 6:
                if(r < 2) return -(254268801024L - 93764321280L*Power(r,2) + 16504750080L*Power(r,4) - 
                        1880186880*Power(r,6) + 187799040*Power(r,8) - 28828800*Power(r,9) - 
                        2515968*Power(r,10) + 655200*Power(r,11) + 6720*Power(r,12) - 3600*Power(r,13) + 
                        5*Power(r,15))/(4.42368e6*(-270257 + 905418*ArcCoth(5) + 6245*Math.log(2)));
                else if(r < 4) return (55251763200L - 768333447168L*r + 539711078400L*Power(r,2) - 458417111040L*Power(r,3) + 
                        484016332800L*Power(r,4) - 260453007360L*Power(r,5) + 59779399680L*Power(r,6) + 
                        509583360*Power(r,7) - 3113510400L*Power(r,8) + 563397120*Power(r,9) - 
                        9609600*Power(r,10) - 7547904*Power(r,11) + 655200*Power(r,12) + 
                        20160*Power(r,13) - 3600*Power(r,14) + 5*Power(r,16))/
                      (8.84736e6*r*(-270257 + 905418*ArcCoth(5) + 6245*Math.log(2)));
                else return -(Power(-6 + r,10)*(-66240 + 5184*r + 35280*Power(r,2) + 11040*Power(r,3) + 
                        1260*Power(r,4) + 60*Power(r,5) + Power(r,6)))/
                        (8.84736e6*r*(-270257 + 905418*ArcCoth(5) + 6245*Math.log(2)));
            default: throw new IllegalArgumentException();
        }
    }
    
    /**
     * Gaussian approximation to value returned by separationProbability.
     */
    public double separationProbabilityGaussian(int n, double r) {

        double b;

        switch(n) {
            case 2: b = 1.19776;
                break;
            case 3: b = 0.674043;
                break;
            case 4: b = 0.543464;
                break;
            case 5: b = 0.4489;
                break;
            case 6: b = 0.381788;
                break;
            case 7: b = 0.331708;
                break;
            case 8: b = 0.29308;
            default: throw new IllegalArgumentException();
        }
        
        return (2*Math.sqrt(b/Math.PI)) * Math.exp(-b*r*r);

    }

    
    private double Power(double x, int n) {
        switch(n) {
        case 2: return x*x;
        case 3: return x*x*x;
        case 4: 
            double x2 = x*x;
            return x2*x2;
        case 5:
            x2 = x*x;
            return x2*x2*x;
        case 6:
            double x3 = x*x*x;
            return x3*x3;
        case 7:
            x3 = x*x*x;
            return x3*x3*x;
        case 8:
            double x4 = x*x*x*x;
            return x4*x4;
        case 9:
            x4 = x*x*x*x;
            return x4*x4*x;
        case 10:
            x2 = x*x;
            x3 = x2*x;
            double x5 = x3*x2;
            return x5*x5;
        default:
            double xp = x;
            for(int j=1; j<n; j++) {
                xp *= x;
            }
            return xp;
        }
    }
    
    private double ArcCoth(double x) {
        return 0.5 * Math.log((x+1)/(x-1));
    }

    public static void main(String[] args) {
        
        for(int n=2; n<10; n++) {
//            ClusterChainWheatley cc = new ClusterChainWheatley(n, null);
            ClusterSinglyConnected cs = new ClusterSinglyConnected(n, null);
            ClusterChain cr = new ClusterChain(n, null, true);
            ClusterChain cc2 = new ClusterChain(n, null);
            System.out.println(n+"\t"+cs.numDiagrams()+
                    "\t"+cr.numDiagrams()+"\t"+cc2.numDiagrams());
            
        }
    }
}
