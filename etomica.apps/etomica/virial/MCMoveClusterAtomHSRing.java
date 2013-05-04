package etomica.virial;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;
import etomica.space3d.RotationTensor3D;
import etomica.util.RandomMersenneTwister;
import etomica.util.RandomNumberGeneratorUnix;

public class MCMoveClusterAtomHSRing extends MCMoveAtom {

    protected final double biasdr;
    public MCMoveClusterAtomHSRing(IRandom random, ISpace _space, double sigma) {
        super(random, null, _space);
        this.sigma = sigma;
        rotAxis = space.makeVector();
        axis0 = space.makeVector();
        axis1 = space.makeVector();
        axis2 = space.makeVector();
        rotTensor = new RotationTensor3D();
        double[] b = new double[]{0,1.1,0.9*1.19776,0.674043,0.543464,0.4489,0.381788,0.331708};
        gaussPrefac = new double[b.length];
        bias = new double[b.length*2+1][0];
        normalWidth = new double[b.length];
        numInserts = new long[b.length*2+1];
        numTrials = new long[b.length*2+1];
        biasdr = 0.1;
        for (int i=1; i<b.length; i++) {
            normalWidth[i] = Math.sqrt(0.5/b[i]);
            
            // (2*Math.sqrt(b/Math.PI))
            gaussPrefac[i] = 2*Math.sqrt(b[i]/Math.PI);
        }
        for (int i=1; i<bias.length; i++) {
            bias[i] = new double[(int)Math.round(i/biasdr)+1];
        }
        twoOSqrtPI = 2/Math.sqrt(Math.PI);
        twoOSqrtPI3 = twoOSqrtPI*twoOSqrtPI*twoOSqrtPI;
    }
    
    public void setBox(IBox box) {
        super.setBox(box);
        inserted = new boolean[box.getLeafList().getAtomCount()];
    }

    public boolean doTrial() {
        
        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.getAtomCount();

        leafAtoms.getAtom(0).getPosition();

        inserted[0] = true;
        for (int i=1; i<n; i++) {
            inserted[i] = false;
        }

        totalCount++;
        while (true) {
            int prevInserted = 0;
            boolean didInsert = false;
            for (int i=1; i<n+1; i++) {
                if (i==n || inserted[i]) {
                    if (i>prevInserted+1) {
                        numInserts[i-prevInserted]++;
                        int nextInserted = i;
                        if (i==n) nextInserted=0;
                        int j = (prevInserted + i)/2;
                        int xj = j-prevInserted;
                        if (xj>2 && (xj&1)==1 && xj==(i-j)) {
                            // we divided an even gap into two odd ones
                            j--;
                        }
                        inserted[j] = true;
                        IVectorMutable pos = leafAtoms.getAtom(j).getPosition();
                        // insert j (between prevInserted and nextInserted)
                        if (i-prevInserted == 2) {
                            // insertion into a lens
                            IVector posPrev = leafAtoms.getAtom(prevInserted).getPosition();
                            IVectorMutable posNext = leafAtoms.getAtom(nextInserted).getPosition();
                            
                            axis0.Ev1Mv2(posNext, posPrev);
                            double r2 = axis0.squared();
                            double rij = Math.sqrt(r2);
                            axis0.TE(1.0/rij);

                            if (Math.abs(axis0.getX(0)) > 0.9999999) {
                                axis1.setX(0,0);
                                axis1.setX(1,1);
                                axis1.setX(2,0);

                                axis2.setX(0,0);
                                axis2.setX(1,0);
                                axis2.setX(2,1);
                            }
                            else if (Math.abs(axis0.getX(0)) > 0.6) {
                                rotAxis.setX(0,1);
                                rotAxis.setX(1,0);
                                rotAxis.setX(2,0);
                                rotAxis.XE(axis0);
                                rotAxis.normalize();
                                rotTensor.setRotationAxisCT(rotAxis, axis0.getX(0));
                                
                                axis1.setX(0,0);
                                axis1.setX(1,1);
                                axis1.setX(2,0);
                                rotTensor.transform(axis1);
                                
                                axis2.setX(0,0);
                                axis2.setX(1,0);
                                axis2.setX(2,1);
                                rotTensor.transform(axis2);
                            }
                            else if (Math.abs(axis0.getX(1)) > 0.6)  {
                                rotAxis.setX(0,0);
                                rotAxis.setX(1,1);
                                rotAxis.setX(2,0);
                                rotAxis.XE(axis0);
                                rotAxis.normalize();
                                rotTensor.setRotationAxisCT(rotAxis, axis0.getX(1));
                                
                                axis1.setX(0,1);
                                axis1.setX(1,0);
                                axis1.setX(2,0);
                                rotTensor.transform(axis1);
                                
                                axis2.setX(0,0);
                                axis2.setX(1,0);
                                axis2.setX(2,1);
                                rotTensor.transform(axis2);
                            }
                            else {
                                rotAxis.setX(0,0);
                                rotAxis.setX(1,0);
                                rotAxis.setX(2,1);
                                rotAxis.XE(axis0);
                                rotAxis.normalize();
                                rotTensor.setRotationAxisCT(rotAxis, axis0.getX(2));
                                
                                axis1.setX(0,1);
                                axis1.setX(1,0);
                                axis1.setX(2,0);
                                rotTensor.transform(axis1);
                                
                                axis2.setX(0,0);
                                axis2.setX(1,1);
                                axis2.setX(2,0);
                                rotTensor.transform(axis2);
                            }

                            double max0 = 2 - rij;
                            double max12 = Math.sqrt(4 - r2);
                            while (true) {
                                numTrials[i-prevInserted]++;
                                pos.Ev1Pv2(posPrev, posNext);
                                pos.TE(0.5);
                                pos.PEa1Tv1(max0*(random.nextFixedDouble()-0.5), axis0);
                                pos.PEa1Tv1(max12*(random.nextFixedDouble()-0.5), axis1);
                                pos.PEa1Tv1(max12*(random.nextFixedDouble()-0.5), axis2);
                                double r2p = pos.Mv1Squared(posPrev);
                                if (r2p > 1) continue;
                                double r2n = pos.Mv1Squared(posNext);
                                if (r2n < 1) break;
                            }
                        }
                        else if (false && i-prevInserted == 3) {
                            double sn = normalWidth[2];
                            IVector posPrev = leafAtoms.getAtom(prevInserted).getPosition();
                            IVector posNext = leafAtoms.getAtom(nextInserted).getPosition();
                            while (true) {
                                numTrials[i-prevInserted]++;
                                double x = random.nextGaussian();
                                double x2 = x*x;
                                pos.setX(0, x*sn);
                                x = random.nextGaussian();
                                x2 += x*x;
                                pos.setX(1, x*sn);
                                x = random.nextGaussian();
                                x2 += x*x;
                                pos.setX(2, x*sn);
                                pos.PE(posNext);
                                if (pos.Mv1Squared(posPrev) < 1) {
                                    // compare Gaussian with true probability
                                    double pGauss = twoOSqrtPI3* Math.exp(-0.5*x2);
                                    double rn = Math.sqrt(pos.Mv1Squared(posNext));
                                    double pn = separationProbability(i-j,rn);
                                    if (pn == 0) continue;
                                    if (pn/pGauss>bias[i-prevInserted][0]) {
                                        bias[i-prevInserted][0] = pn/pGauss;
                                    }
                                    if (random.nextDouble() < pn/(bias[i-prevInserted][0]*pGauss)) {
                                        break;
                                    }                    
                                }
                            }
                        }
                        else {
                            // Gaussian distributions with bisection
                            IVector posPrev = leafAtoms.getAtom(prevInserted).getPosition();
                            IVector posNext = leafAtoms.getAtom(nextInserted).getPosition();
                            double rpn = Math.sqrt(posNext.Mv1Squared(posPrev));
                            while (true) {
                                numTrials[i-prevInserted]++;
                                double sp = normalWidth[j-prevInserted];
                                double sp2 = sp*sp;
                                double sn = normalWidth[i-j];
                                double sn2 = sn*sn;
                                double s = Math.sqrt(sp2*sn2/(sp2+sn2));
                                double x = random.nextGaussian();
                                double x2 = x*x;
                                pos.setX(0, x*s);
                                x = random.nextGaussian();
                                x2 += x*x;
                                pos.setX(1, x*s);
                                x = random.nextGaussian();
                                x2 += x*x;
                                pos.setX(2, x*s);
                                // (2*Math.sqrt(b/Math.PI))
                                // (2/(sqrt(b/pi))^3
                                double pGauss = gaussPrefac[j-prevInserted]*gaussPrefac[i-j] * Math.exp(-0.5*x2);
                                pos.PEa1Tv1(sn2/(sp2+sn2), posPrev);
                                pos.PEa1Tv1(sp2/(sp2+sn2), posNext);
                                double rp = Math.sqrt(pos.Mv1Squared(posPrev));
                                double rn = Math.sqrt(pos.Mv1Squared(posNext));
                                double pp = separationProbability(j-prevInserted,rp);
                                if (pp == 0) continue;
                                double pn = separationProbability(i-j,rn);
                                if (pn == 0) continue;
                                double p = pp*pn;
                                int biasBin = (int)(rpn/biasdr);
                                if (p/pGauss>bias[i-prevInserted][biasBin]) {
//                                    if (i-prevInserted==8) {
//                                        System.out.println((i-prevInserted)+" "+rpn+" "+Math.sqrt(x2)+" "+p/pGauss);
  //                                  }
                                    bias[i-prevInserted][biasBin] = p/pGauss;
                                }
                                if (random.nextDouble() < p/(bias[i-prevInserted][biasBin]*pGauss)) {
                                    break;
                                }
                            }
                        }
                        didInsert = true;
                    }
                    prevInserted = i;
                }
            }
            if (!didInsert) break;
        }
//        for (int i=1; i<n; i++) {
//            IVector pos = leafAtoms.getAtom(i).getPosition();
//            System.out.println(Math.sqrt(pos.Mv1Squared(leafAtoms.getAtom(i-1).getPosition())));
//        }
//        System.out.println(Math.sqrt(leafAtoms.getAtom(0).getPosition().Mv1Squared(leafAtoms.getAtom(n-1).getPosition())));

        if (false && random.nextDouble() < 0.0001) {
            long s = 0;
            for (int i=0; i<numTrials.length; i++){
                if (numInserts[i] > 0) {
                    s += numTrials[i];
                    System.out.print(String.format("%2d  %5.3f  %5.3f\n", i, ((double)numInserts[i])/numTrials[i], ((double)numTrials[i])/totalCount));
//                    boolean something = false;
//                    for (int j=0; j<bias[i].length; j++) {
//                        if (bias[i][j] == 0) {
//                            if (something) break;
//                            continue;
//                        }
//                        System.out.print(String.format("  %4.2f  %6.4f\n", j*biasdr, bias[i][j]));
//                        something = true;
//                    }
                }
            }
            System.out.println(" total "+(s/(((double)n)*totalCount)));
            System.out.println();
        }
		((BoxCluster)box).trialNotify();
		return true;
	}
	
    /**
     * Normalized probability for the separation of the centers of two spheres located at the
     * end of a chain of (n+1) overlapping spheres. All spheres of unit diameter.
     * @param n number of bonds in the chain, such that number of spheres is n+1
     * @param r end-to-end separation
     * @return value of normalized probability density
     */
    public static double separationProbability(int n, double r) {

        if(r >= n) return 0;

        switch(n) {
            case 1:
                return 1;
            case 2:
                double rm2 = r-2;
                return rm2*rm2*(4+r)/12;
            case 3:
                if(r < 1) {
                    double r2 = r*r;
                    return (-525 + (315 - (63 - r2)*r2)*r2)/fac3;
                }
                double rm3 = r-3;
                double rm32 = rm3*rm3;
                return -rm32*rm32*(-6 + (27 + (12 + r)*r)*r)/(2*r*fac3);
            case 4:
                if(r < 2) {
                    double r2 = r*r;
                    return (-348160 + r2*(184320 + r2*(-56448 + r*(15120 + 960*r + 
                        r2*(-540 + 3*r2)))))/fac4;
                }
                double rm4 = r-4;
                double rm43 = rm4*rm4*rm4;
                return -(rm43*rm43*(-144 + r*(224 + r*(156 + r*(24 + r)))))/(fac4*r);
            case 5:
                if(r < 1) {
                    double r2 = r*r;
                    return -(146392675 + r2*(-63050130 + r2*(12657645 + r2*(-1501500 + 
                        r2*(96525 + r2*(-1170 + r2*3))))))/(2*fac5);
                }
                if(r < 3) return (335430 + r*(-75822175+ r*(8925150 + r*(14401530 + 
                        r*(20091500 + r*(-20765745 + r*(5791500 + 
                        r*(-42900 + r*(-244530 + r*(32175 + r*(1430 + 
                        r*(-390 + r*r))))))))))))/(r*fac5);
                double rm5 = r-5;
                double rm52 = rm5*rm5;
                double rm54 = rm52*rm52;
                return -(rm54*rm54*(-3060 + r*(1825 + r*(2260 + r*(510 + r*(40 + r))))))/(4*r*fac5);
            case 6:
                if(r < 2) {
                    double r2 = r*r;
                    return -(254268801024.0 + r2*(-93764321280.0 + r2*(16504750080.0 + 
                        r2*(-1880186880 + r2*(187799040 + r*(-28828800 + 
                        r*(-2515968 + r*(655200 + r*(6720 + r*(-3600 + 
                        r2*5))))))))))/fac6;
                }
                if(r < 4) return (55251763200.0 + r*(-768333447168.0 + r*(539711078400.0 + r*(-458417111040.0 + 
                        r*(484016332800.0 + r*(-260453007360.0 + r*(59779399680.0 + 
                        r*(509583360 + r*(-3113510400.0 + r*(563397120 +
                        r*(-9609600 + r*(-7547904 + r*(655200 + r*(20160 + r*(-3600 + r*r*5)))))))))))))))/(2*r*fac6);
                double rm6 = r-6;
                double rm62 = rm6*rm6;
                double rm64 = rm62*rm62;
                return -(rm64*rm64*rm62*(-66240 + r*(5184 + r*(35280 + r*(11040 + 
                        r*(1260 + r*(60 + r)))))))/(2*r*fac6);
            case 7:
                if(r < 1) {
                    double r2 = r*r;
                    return -2.5163301439885776e12 + 
                        r2*(809442607918.0 + r2*(-124848195160.0 + 
                         r2*(1.2293836678933332e10 + 
                            r2*(-8.612955208e8 + 
                               r2*(4.37602984e7 + 
                                  r2*(-1.4048346666666667e6 + 
                                     r2*(12920 + r2*(-38 + (2*r2)/63.))))))))/fac7;
                }
                if(r < 3) {
                    return (8068545450.0 + r*(-1585383868067349.0 + r*(521067793680.0 + 
                            r*(508237500736665.0 + r*(3778612050360.0 + r*(-84552506356500.0 + r*(6652988339760.0 + 
                               r*(2312456534388.0 + r*(3149656883100.0 + r*(-1770858823734.0 + r*(277928450800.0 + 
                               r*(12032729982.0 + r*(-8010180360.0 + r*(572620860 + r*(47209680 + 
                               r*(-6104700 + r*(-67830 + r*(17955 - 15*r*r))))))))))))))))))/(630.*r*fac7);
                }
                if(r < 5) {
                    return (562806119671740.0 + r*(-2815813630010031.0 + r*(3237357537947232.0 + 
                            r*(-2748586830655365.0 + r*(1748138067798864.0 + 
                            r*(-667883097963900.0 + r*(113164463541024.0 + 
                            r*(11351302774812.0 + r*(-9238129369560.0 + 
                            r*(1660864339134.0 + r*(-59025846880.0 + 
                            r*(-22032228582.0 + r*(3204072144.0 + 
                            r*(-59826060 + r*(-18883872 + r*(1220940 + r*(27132 + r*(-3591 + 3*r*r)))))))))))))
                                    )))))/(315.*r*fac7);
                }
                return (20675145930749730.0 + r*(-29368263246348187.0 + r*(9700068982144464.0 + 
                        r*(4733961739079895.0 + r*(-4668230192527272.0 + 
                        r*(1459645883759700.0 + r*(-153103614167952.0 + 
                        r*(-28887620637876.0 + r*(11267620198380.0 + 
                        r*(-1279561765482.0 + r*(-13941687760.0 + 
                        r*(18247233186.0 + r*(-1602036072 + 
                        r*(-10445820 + r*(9441936 + r*(-406980 + r*(-13566 + r*(1197 - r*r)))))))))))))))))
                        )/(630.*r*fac7);
//                is this preferable?
//                double rm7 = (r-7)*(r-7);
//                rm7 = rm7*rm7*rm7;//(r-7)^6  need (r-7)^12
//                return -rm7*rm7*(-1493730 + r*(-438893 + r*(558768 + r*(248955 + r*(37870 + r*(2625 + r*(84 + r)))))))/(630.*r*fac7);
            default: throw new IllegalArgumentException(); 
        }
    }

    private static double Power(double x, int n) {
        double xp = x;
        for(int j=1; j<n; j++) {
            xp *= x;
        }
        return xp;
    }
    
    private static double ArcCoth(double x) {
        return 0.5 * Math.log((x+1)/(x-1));
    }

    public double getA() {
        return 1;
    }

    public double getB() {
    	return 0.0;
    }
    
    public void rejectNotify() {
        throw new RuntimeException("nope");
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    protected final double sigma;
    protected final RotationTensor3D rotTensor;
    protected final IVectorMutable rotAxis, axis0, axis1, axis2;
    protected final double[] normalWidth;
    protected final double[][] bias;
    protected final double twoOSqrtPI, twoOSqrtPI3;
    protected final double[] gaussPrefac;
    protected boolean[] inserted;
    protected long[] numInserts, numTrials;
    protected long totalCount;
    protected final static double fac3 = (9.*(-92 + 27*Math.log(3.)));
    protected final static double fac4 = (4608.*(-179 + 256*ArcCoth(3)));
    protected final static double fac5 = (90.*(-2773712 + 6640625*ArcCoth(4) + 3727*Math.log(3)));
    protected final static double fac6 = (4.42368e6*(-270257 + 905418*ArcCoth(5) + 6245*Math.log(2)));
    protected final static double fac7 = -15825661837824. + 3573372188392.*ArcCoth(4) + 65635383907142.*ArcCoth(6) + 12807215*Math.log(3);
    
    public static void main1(String[] args) {
        if (false) {
            double dr = 0.01;
            for (int i=2; i<7; i++) {
                for (int jx=1; jx<100*i; jx++) {
                    double r= dr*jx;
                    System.out.println(r+" "+separationProbability(i, r));
                }
                System.out.println("&");
            }
        }
        double[] b = new double[]{0,0.4,0.9*1.19776,0.674043,0.543464,0.4489,0.381788,0.331708};
        RandomMersenneTwister random = new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray());
        FileWriter fw = null;
//        try {
//            fw = new FileWriter("maxbias.dat");
//        }
//        catch (IOException e) {
//            throw new RuntimeException(e);
//        }
        for (int i=3; i<=12; i++) {
            int j = i/2;
            double normalWidthP = Math.sqrt(0.5/b[j]);
            double normalWidthN = Math.sqrt(0.5/b[i-j]);
            double gaussPrefacP = 2*Math.sqrt(b[j]/Math.PI);
            double gaussPrefacN = 2*Math.sqrt(b[i-j]/Math.PI);
            double bias = 0;
            int lastBiasBin = -1;
            for (int jr=0; jr<100*i; jr++) {
                double rij = jr*0.01;
                int biasBin = jr/10;
                if (biasBin > lastBiasBin) {
                    if (lastBiasBin > -1) {
                        System.out.println(i+" "+lastBiasBin*0.1+" "+bias);
                        if (fw != null) {
                            try {
                                fw.write(i+" "+lastBiasBin+" "+bias+"\n");
                            }
                            catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        }
                    }
                    bias = 0;
                }
                double sp = normalWidthP;
                double sp2 = sp*sp;
                double sn = normalWidthN;
                double sn2 = sn*sn;
                double s = Math.sqrt(sp2*sn2/(sp2+sn2));
                double cx = sp2/(sp2+sn2)*rij;
                long nk = 100000;
                for (long k=0; k<=nk; k++) {
                    double x = random.nextGaussian();
                    if (k==0 || (i&1)==0) x = 0;
                    double x2 = x*x;
                    double jx = x*s;
                    x = random.nextGaussian();
                    if (k==0) x = 0;
                    x2 += x*x;
                    double jy = x*s;
                    x = random.nextGaussian();
                    if (k==0) x = 0;
                    x2 += x*x;
                    double jz = x*s;
                    // (2*Math.sqrt(b/Math.PI))
                    // (2/(sqrt(b/pi))^3
                    double pGauss = gaussPrefacP*gaussPrefacN * Math.exp(-0.5*x2);
                    jx += cx;
                    double rp = Math.sqrt(jx*jx + jy*jy + jz*jz);
                    double rn = Math.sqrt((rij-jx)*(rij-jx) + jy*jy + jz*jz);
                    double pp = separationProbability(j,rp);
                    if (pp == 0) continue;
                    double pn = separationProbability(i-j,rn);
                    if (pn == 0) continue;
                    double p = pp*pn;
                    if (p/pGauss>bias) {
                        System.out.println(i+" "+rij+" "+(jx-cx)+" "+Math.sqrt(jy*jy+jz*jz)+" "+p/pGauss);
                        bias = p/pGauss;
                    }
                }
                lastBiasBin = biasBin;
            }
            if (fw != null) {
                System.out.println(i+" "+lastBiasBin*0.1+" "+bias);
                try {
                    fw.write(i+" "+lastBiasBin+" "+bias+"\n");
                }
                catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
        if (fw != null) {
            try {
                fw.close();
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
