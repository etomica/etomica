/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;
import etomica.virial.BoxCluster;

import java.io.FileWriter;
import java.io.IOException;

public class MCMoveClusterAtomHSRing extends MCMoveBox {

    protected final IRandom random;
    private static final int nLookUp = 5000;
    private static final double[] lookUpL;
    protected final Vector3D standardLensPoint, normalVector;
    protected double[] maxBias, maxBiasR, p0;
    
    public int[] seq;
    
    static {
        lookUpL = new double[nLookUp];
        for(int i=0; i<nLookUp; i++) {
            double U = (double)(i+1)/nLookUp/1.5;
            lookUpL[i] = lensRootExact(U);
        }
    }
    
    private static final double lensRootExact(double U) {
        double theta = (Math.atan2((Math.sqrt(3)*Math.sqrt((4 - 3*U)*U)),(-2 + 3*U)))/3.;
        return 1 + Math.cos(theta) - Math.sqrt(3)*Math.sin(theta);
    }
    
    private static final double lensRootApprox(double U) {
        double sqrtU = Math.sqrt(U);
        return sqrtU*(1 + sqrtU*(0.16666666666666666 + 
                sqrtU*(0.06944444444444445 + 
                        sqrtU*(0.037037037037037035 + 
                           sqrtU*(0.02228009259259259 + 
                              sqrtU*(0.01440329218106996 + 
                                 (0.009769643775720165 + (5*sqrtU)/729.)*sqrtU))))));
    }

    public MCMoveClusterAtomHSRing(IRandom random, Box box, double sigma) {
        super();
        this.random = random;
        this.sigma = sigma;
        setBox(box);
        axis0 = box.getSpace().makeVector();
        standardLensPoint = (Vector3D)box.getSpace().makeVector();
        normalVector = (Vector3D)box.getSpace().makeVector();
        double[] b = new double[]{0,1.1,0.9*1.19776,0.674043,0.543464,0.4489,0.381788,0.331708,0.29308,0.262442,0.237574};
        normalWidth = new double[b.length];
        numInserts = new long[b.length*2+1];
        numTrials = new long[b.length*2+1];
        for (int i=1; i<b.length; i++) {
            normalWidth[i] = Math.sqrt(0.5/b[i]);
        }
        maxBias = new double[b.length*2+1];
        maxBiasR = new double[b.length*2+1];
        p0 = new double[b.length*2+1];
        // this approach to computing maxBias only works for symmetric insertions
        // in doing the trial, we avoid asymmetric insertions with more than one
        // point on each side by doing (for instance) 1+6 (handled without using
        // maxBias computed here) instead of 3+4 (which would need the maxBias
        // computed here).
        for (int i=2; i<b.length; i++) {
            p0[i] = separationProbability(i, 0);
            double sp2 = 0.5/b[i];
            double s = 0.5/Math.sqrt(b[i]);
            double rmin = 0;
            double maxr = -1;
            double maxbias = -1;
            for (int k=5; k<=50; k+=3) {
                double dr = s/(1L<<k);
                for (int j=0; j<1000; j++) {
                    double r = rmin + j*dr;
                    double pAct = separationProbability(i, r);
                    pAct *= pAct;
                    double pGauss = Math.exp(-r*r/sp2);
                    double pRatio = pAct/pGauss;
                    if (pRatio > maxbias) {
                        maxbias = pRatio;
                        maxr = r;
                    }
                }
//                System.out.print(String.format("%2d %2d %22.15e %19.15f %19.15f\n", i, k, dr, maxr, maxbias));
                rmin = maxr - 500*dr/8;
            }
            maxBias[i] = maxbias;
            maxBiasR[i] = maxr;
        }
    }

    public void setBox(Box box) {
        super.setBox(box);
        inserted = new boolean[box.getLeafList().size()];
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
    }

    public boolean doTrial() {
        
        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.size();

        leafAtoms.get(0).getPosition();

        inserted[0] = true;
        for (int i=1; i<n; i++) {
            inserted[i] = false;
        }

        totalCount++;

        if (seq == null) {
            seq = new int[n];
        }
        for (int i=0; i<n; i++) {
            seq[i] = i;
        }
        for (int i=0; i<n; i++) {
            int j = i+random.nextInt(n-i);
            int k = seq[j];
            seq[j] = seq[i];
            seq[i] = k;
        }
//        System.out.println(Arrays.toString(seq));
        leafAtoms.get(seq[0]).getPosition().E(0);
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
                        if (((i-prevInserted)&1)==1) {
                            j = prevInserted+1;
                        }
                        inserted[j] = true;
                        Vector pos = leafAtoms.get(seq[j]).getPosition();
                        // insert j (between prevInserted and nextInserted)
                        if (i-prevInserted == 2) {

                            // insertion into a lens
                            Vector posPrev = leafAtoms.get(seq[prevInserted]).getPosition();
                            Vector posNext = leafAtoms.get(seq[nextInserted]).getPosition();
                            
                            axis0.Ev1Mv2(posNext, posPrev);
                            double r2 = axis0.squared();
                            if (r2 < 1) {
                                randomLensPointInBox(posPrev, posNext, pos);
                            }
                            else {
                                randomLensPoint(posPrev, posNext, (Vector3D)pos, axis0, r2);
                            }
                        }
                        else if (j-prevInserted==1) {
                            // odd gap.  even it up  n = 1 + (n-1)
                            Vector posPrev = leafAtoms.get(seq[prevInserted]).getPosition();
                            Vector posNext = leafAtoms.get(seq[nextInserted]).getPosition();
                            while (true) {
                                numTrials[i-prevInserted]++;
                                pos.setRandomInSphere(random);
                                pos.PE(posPrev);
                                double rn = Math.sqrt(pos.Mv1Squared(posNext));
                                double pn = separationProbability(i-j,rn);
                                if (pn == 0) continue;
                                double p = pn;
                                double rand = random.nextDouble();

                                double r2 = posPrev.Mv1Squared(posNext);
                                double fullBias = 0;
                                if (r2 < 1) {
                                    // if prev core overlaps next center, then maxBias is the separationProbability(0)
                                    fullBias = p0[i-j];
                                }
                                else {
                                    // maxBias is separationProbability evaluated at the edge of prev's core
                                    double myR = Math.sqrt(r2) - 1;
                                    fullBias = separationProbability(i-j, myR);
                                }
                                if (rand < p/fullBias) {
                                    break;
                                }
                            }
                        }
                        else {
                            // Gaussian distributions with bisection
                            if (j-prevInserted != i-j) {
                                throw new RuntimeException("maxBias has only been computed for symmetric insertion");
                            }
                            Vector posPrev = leafAtoms.get(seq[prevInserted]).getPosition();
                            Vector posNext = leafAtoms.get(seq[nextInserted]).getPosition();
                            double rpn = Math.sqrt(posNext.Mv1Squared(posPrev));
                            while (true) {
                                numTrials[i-prevInserted]++;
                                double sp = normalWidth[j-prevInserted];
                                double sp2 = sp*sp;
                                double sn = normalWidth[i-j];
                                double sn2 = sn*sn;
                                double s = Math.sqrt(sp2*sn2/(sp2+sn2));
                                double x = random.nextGaussian();
                                pos.setX(0, x*s);
                                x = random.nextGaussian();
                                pos.setX(1, x*s);
                                x = random.nextGaussian();
                                pos.setX(2, x*s);
                                // (2*Math.sqrt(b/Math.PI))
                                // (2/(sqrt(b/pi))^3
                                pos.PEa1Tv1(sn2/(sp2+sn2), posPrev);
                                pos.PEa1Tv1(sp2/(sp2+sn2), posNext);
                                double rp = Math.sqrt(pos.Mv1Squared(posPrev));
                                double rn = Math.sqrt(pos.Mv1Squared(posNext));
                                double pp = separationProbability(j-prevInserted,rp);
                                if (pp == 0) continue;
                                double pn = separationProbability(i-j,rn);
                                if (pn == 0) continue;
                                double p = pp*pn;

                                double fullBias = 0;
                                if (rpn < 2*maxBiasR[i-j]) {
                                    // if prev and next are close enough, the maxBias is our precomputed value
                                    fullBias = maxBias[i-j];
                                }
                                else {
                                    // maxBias is p ratio evaluated at the midpoint
                                    fullBias = separationProbability(i-j, 0.5*rpn);
                                    fullBias *= fullBias;
                                    fullBias /= Math.exp(-0.5*(0.5*rpn*0.5*rpn + 0.5*rpn*0.5*rpn)/sp2);
                                }
                                // recalc without normalization
                                double pGaussFull = Math.exp(-0.5*(rp*rp/sp2 + rn*rn/sn2));
                                if (random.nextDouble() < p/(fullBias*pGaussFull)) {
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
        if (sigma != 1) {
            for (int i=1; i<n; i++) {
                leafAtoms.get(seq[i]).getPosition().TE(sigma);
            }
        }
//        System.out.println(Math.sqrt(leafAtoms.getAtom(0).getPosition().Mv1Squared(leafAtoms.getAtom(n-1).getPosition())));

        if (false && random.nextDouble() < 0.00001) {
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
     * Sets a given vector to a point chosen uniformly at random in the lens-shaped region
     * formed by the overlap of two spheres of unit radius located at the given positions r1 and r2.
     * @param r1 (input) position of one sphere
     * @param r2 (input) position of other sphere
     * @param point (output) random position selected in overlap of the given spheres
     * @throws IllegalArgumentException if given spheres do not overlap
     */
    protected void randomLensPoint(Vector r1, Vector r2, Vector3D point, Vector r12, double rSquared) {
        double d = Math.sqrt(rSquared); //distance between spheres
        randomLensPoint(d, standardLensPoint); //select point for two spheres in standard configuration and separated by this amount
        r12.TE(1./d);//normalize
        normalVector.setPerpendicularTo(r12);
        normalVector.normalize();
        point.Ev1Pv2(r1,r2);
        point.TE(0.5);//midpoint between spheres
        point.PEa1Tv1(standardLensPoint.getX(0),r12);
        point.PEa1Tv1(standardLensPoint.getX(1),normalVector);
        normalVector.XE(r12);
        point.PEa1Tv1(standardLensPoint.getX(2),normalVector);
    }

    protected void randomLensPointInBox(Vector r0, Vector r1, Vector point) {
        double minX = r0.getX(0) - 1;
        double a = r1.getX(0) - 1;
        double maxX = 0;
        if (a>minX) {
            maxX = minX+2;
            minX = a;
        }
        else {
            maxX = a+2;
        }
        double minY = r0.getX(1) - 1;
        a = r1.getX(1) - 1;
        double maxY = 0;
        if (a>minY) {
            maxY = minY+2;
            minY = a;
        }
        else {
            maxY = a+2;
        }
        double minZ = r0.getX(2) - 1;
        a = r1.getX(2) - 1;
        double maxZ = 0;
        if (a>minZ) {
            maxZ = minZ+2;
            minZ = a;
        }
        else {
            maxZ = a+2;
        }
        double dx = maxX - minX;
        double dy = maxY - minY;
        double dz = maxZ - minZ;
        while (true) {
            point.setX(0, minX + dx*random.nextFixedDouble());
            point.setX(1, minY + dy*random.nextFixedDouble());
            point.setX(2, minZ + dz*random.nextFixedDouble());
            if (point.Mv1Squared(r0) < 1 && point.Mv1Squared(r1) < 1) return;
        }            
    }

    /** 
     * Returns a randomly selected point from within the lens-shaped region of
     * overlap of two spheres of unit radius separated by a distance d.  The positions
     * of the spheres are on the x-axis, with the midpoint between them at the origin.
     * 
     * This is used to select a point where a third sphere could be positioned such that
     * it would overlap the two others, with all having unit diameter.
     * 
     * @param d the distance between the spheres
     * @param point (output) the randomly selected point
     * @throws IllegalArgumentException if d > 2, corresponding to non-overlapping spheres
     */
    protected void randomLensPoint(double d, Vector3D point) {
        
        if(d > 2) throw new IllegalArgumentException("No overlap point if spheres are separated by a distance "+d);
        
        double H = 1.0 - 0.5*d; //distance from origin to tip of sphere
        double U = random.nextFixedDouble() * H*H*(1.-H/3.);
        double h = lensRoot(U); //randomly-selected distance from tip of sphere 
        double x = H - h; //distance from center plane of sphere pair (origin)
        double xS = 1 - h; //distance from center of negative-x sphere
        if(random.nextFixedDouble() > 0.5) x = -x;
        double circleDiameter = 2*Math.sqrt(1 - xS*xS);
        
        double y, z;
        do {
            y = random.nextFixedDouble() - 0.5;
            z = random.nextFixedDouble() - 0.5;
        } while(y*y + z*z > 0.25);
        
        y *= circleDiameter;
        z *= circleDiameter;

        point.E(x, y, z);
        
    }
    
    // desired root of a^2 - a^3/3 - U = 0, via look-up and Newton iteration
    private static double lensRoot(double U) {
        final double twoThirds = 2./3.; 
        int iR = (int)(Math.round(1.5*U*nLookUp)-1);
        if(iR < 6) return lensRootApprox(U);
        if(iR == nLookUp) iR--;
        
        double a = lookUpL[iR];
        for(int i=0; i<3; i++) {
            //a = (3*a*a - 2*a*a*a + 3*U)/(6*a - 3*a*a);
            a = (a*a*(1 - twoThirds*a) + U)/((2 - a)*a);
        }
        return a;
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
        double rm3, rm32, rm4, rm42, rm5, rm52, rm6, rm62, rm8, rm82;
        double r2 = r*r;

        switch(n) {
            case 1:
                return 1;
            case 2:
                double rm2 = r-2;
                return rm2*rm2*(4+r)/12;
            case 3:
                if(r < 1) {
                    return (-525 + (315 - (63 - r2)*r2)*r2)/fac3;
                }
                rm3 = r-3;
                rm32 = rm3*rm3;
                return -rm32*rm32*(-6 + (27 + (12 + r)*r)*r)/(2*r*fac3);
            case 4:
                if(r < 2) {
                    return (-348160 + r2*(184320 + r2*(-56448 + r*(15120 + 960*r + 
                        r2*(-540 + 3*r2)))))/fac4;
                }
                rm4 = r-4;
                double rm43 = rm4*rm4*rm4;
                return -(rm43*rm43*(-144 + r*(224 + r*(156 + r*(24 + r)))))/(fac4*r);
            case 5:
                if(r < 1) {
                    return -(146392675 + r2*(-63050130 + r2*(12657645 + r2*(-1501500 + 
                        r2*(96525 + r2*(-1170 + r2*3))))))/(2*fac5);
                }
                if(r < 3) {
                    rm3 = r - 3;
                    rm32 = rm3*rm3;
                    return (-11.645519426547745417 + rm3)*(2.9955737735903254375 + rm3)*
                        (8.699501424572235734 + rm3)*(13.97637523661338622 + rm3)*
                        (21.342480273228225428 + rm3)*(0.7342909491741547658 - 1.667110360558664731*rm3 + 
                          rm32)*(1.0303848213615477176 - 1.5307636570742344459*rm3 + rm32)*
                        (2.3897772117330571797 - 0.8886883892789973696*rm3 + rm32)*
                        (15.643325581592968058 + 7.7181511254554691437*rm3 + rm32)/(r*fac5);
                }
                // 3 < r < 5
                rm5 = r-5;
                rm52 = rm5*rm5;
                double rm54 = rm52*rm52;
                return -rm54*rm54*(-0.7806915955741175188 + r)*(2.5006342896385876075 + r)*(6.6436016161512016326 + r)*
                (12.038860037830168394 + r)*(19.597595651954159885 + r)/(4*r*fac5);
            case 6:
                if(r < 2) {
                    return 10*(-20.663272399445758108 + r)*(-12.670131734239847386 + r)*
                        (8.5812267463490890752 + r)*(14.483969838752005272 + r)*(22.472894082217685734 + r)*
                        (12.627138033294569079 - 7.0774518544647587239*r + r2)*
                        (12.791198774472720108 - 6.8418369699588631093*r + r2)*
                        (13.72421926822998938 - 6.1151889532396643285*r + r2)*
                        (6.706456116575184131 + 3.6433022025736211298*r + r2)*
                        (4.6778892336765258878 + 4.1864890414564904446*r + r2)/fac6;
                }
                if(r < 4) { 
                    rm4 = r - 4;
                    rm42 = rm4*rm4;
                    return -5*(-14.155802106465130133 + rm4)*(3.9243026987024033091 + rm4)*
                        (9.1490819706794762747 + rm4)*(13.88047597753425173 + rm4)*
                        (19.812004333376225671 + rm4)*(27.940400019553254824 + rm4)*
                        (0.7636120659947073487 - 1.7189254501991924474*rm4 + rm42)*
                        (0.9408278613331036471 - 1.6531001375241852783*rm4 + rm42)*
                        (1.5061280567943069866 - 1.4416077757701893518*rm4 + rm42)*
                        (3.7569261240788530062 - 0.5759377102201528465*rm4 + rm42)*
                        (20.816947688202424179 + 8.8391081803332382487*rm4 + rm42)/(r*fac6);
                }
                // 4 < r < 6
                rm6 = r-6;
                rm62 = rm6*rm6;
                double rm64 = rm62*rm62;
                return rm64*rm64*rm62*(-1.1092620887890839184 + r)*(2.1435339705616727282 + r)*(6.0415816056834182244 + r)*
                (10.830695965418476141 + r)*(16.896976675134874992 + r)*(25.196473871990641832 + r)/(r*fac6);
            case 7:
                if(r < 1) {
                    double r4 = r2*r2;
                    return -10*(-704.83587167279709762 + r2)*(-323.57493353067866407 + r2)*
                        (-132.54434423563879638 + r2)*(109.33047070592459752 - 20.203460280273764197*r2 + 
                                r4)*(125.05528171708708364 - 15.177992581472341505*r2 + r4)*
                              (191.78332637567910347 - 0.663397699139336226*r2 + r4)/fac7;
                }
                if(r < 3) {
                    return 7.5*(-24.499837185351812732 + r)*(-15.844913345579656722 + r)*
                        (-5.0893323857263886304e-6 + r)*(7.877286483990656681 + r)*
                        (13.136217777622784723 + r)*(19.559680644241281934 + r)*(28.211509940477606668 + r)*
                        (21.151685853233686036 - 9.1828852083594505267*r + r2)*
                        (21.171939154852635151 - 9.0498104825698319163*r + r2)*
                        (21.213151779102138681 - 8.6773249347811283059*r + r2)*
                        (22.077019818148285214 - 7.754572190491953564*r + r2)*
                        (6.853422456448555975 + 2.8324263298735301241*r + r2)*
                        (3.3173879467885562645 + 3.3922272602603593643*r + r2)/(r*fac7);
                }
                if(r < 5) {
                    rm5 = r - 5;
                    rm52 = rm5*rm5;
                    return -3*(-16.662749785120940509 + rm5)*(4.7346182598231468938 + rm5)*
                        (9.6676685774179364482 + rm5)*(14.088235762311330519 + rm5)*
                        (19.329584888899402107 + rm5)*(25.844788727614001316 + rm5)*
                        (34.635638755071043307 + rm5)*(0.7866375529384692204 - 1.7543389558235202893*rm5 + 
                          rm52)*(0.9068337796477158114 - 1.7165744412374368829*rm5 + rm52)*
                        (1.2331424606740560288 - 1.6137441143881413453*rm5 + rm52)*
                        (2.1089511465009502349 - 1.335509407147958674*rm5 + rm52)*
                        (5.4502260558993646031 - 0.2430001912868617859*rm5 + rm52)*
                        (27.11365482880402116 + 10.0253819238679988966*rm5 + rm52)/(r*fac7);
                }
                // 5 < r < 7
                double rm7n = (r-7)*(r-7);//(r-7)^2
                rm7n = rm7n*rm7n*rm7n;//(r-7)^6  need (r-7)^12
                return 0.5*rm7n*rm7n*(-1.4534486570398939108 + r)*(1.7929863166426809281 + r)*(5.5384677768768614701 + r)*
                    (9.9694148838300796512 + r)*(15.30793344875122891 + r)*(21.950316601593546361 + r)*
                    (30.894329629345496591 + r)/(r*fac7);
            case 8:
                if(r<2) {
                    return -35*(-30.582471610870566927 + r)*(-21.4908457310048427 + r)*
                        (-14.484147565373872014 + r)*(10.662597265771953598 + r)*
                        (16.408562447082325391 + r)*(23.288589904948409015 + r)*(32.426989719172574766 + r)*
                        (18.439888836196970696 - 8.5599101758760078894*r + r2)*
                        (18.760761988189045438 - 8.3850301470408463685*r + r2)*
                        (19.658126102966866404 - 7.9403579767720178202*r + r2)*
                        (22.679154244994764765 - 6.9412884434867099229*r + r2)*
                        (14.463130049364862635 + 4.5184148495751038968*r + r2)*
                        (9.813570763629201288 + 5.3799012953213332193*r + r2)*
                        (8.377400849347788983 + 5.6989961685531637556*r + r2)/(fac8);
                    }
                if(r<4) return 21*(-28.311040615879795745 + r)*(-19.041628162010081697 + r)*
                    (-0.0010354301873885649752 + r)*(7.2856078858023037761 + r)*
                    (12.168420318340987233 + r)*(17.827697237859167396 + r)*(24.759841604895890736 + r)*
                    (34.007679826964571043 + r)*(31.754303634910252779 - 11.2611575639418972502*r + r2)*
                    (31.722982965456166222 - 11.1771716842614529753*r + r2)*
                    (31.640203208039229256 - 10.9652290087249278884*r + r2)*
                    (31.383714649410244196 - 10.4723942696851196555*r + r2)*
                    (32.018573532669754486 - 9.3686787745505957835*r + r2)*
                    (7.517238032382953375 + 2.0104551189211066296*r + r2)*
                    (2.2414089287174718245 + 2.5386335164572327466*r + r2)/(fac8*r);
                if(r<6) {
                    rm6 = r - 6;
                    rm62 = rm6*rm6;
                    return -7*(-19.167762398755254177 + rm6)*(5.4699074502141796531 + rm6)*
                        (10.221511396855285446 + rm6)*(14.445389836354267471 + rm6)*
                        (19.27349980745107378 + rm6)*(24.993150191494885028 + rm6)*
                        (32.026278513792358258 + rm6)*(41.40543904812347119 + rm6)*
                        (0.805189915394450896 - 1.7805440363273398668*rm6 + rm62)*
                        (0.8928448055378504483 - 1.7565486615630626029*rm6 + rm62)*
                        (1.110713842408393392 - 1.696817723289354033*rm6 + rm62)*
                        (1.5963919045924004093 - 1.5632022394127791114*rm6 + rm62)*
                        (2.8355860775799921747 - 1.2193665640844762326*rm6 + rm62)*
                        (7.466622130539541303 + 0.1041858474042261884*rm6 + rm62)*
                        (34.399377457727163683 + 11.2448795317425190091*rm6 + rm62)/(fac8*r);
                }
                // 6 < r < 8
                rm8 = r - 8;
                double rm814 = rm8*rm8;//(r-8)^2
                rm814 = rm814*rm814*rm814*rm8;//(r-8)^7
                rm814 = rm814*rm814;//(r-8)^14
                return rm814*(6.1906986749463368728 + rm8)*(9.4433989805794742917 + rm8)*
                    (13.087904646691476115 + rm8)*(17.28304422447324228 + rm8)*
                    (22.179670619491430171 + rm8)*(27.998841449141393342 + rm8)*
                    (35.149253935946578148 + rm8)*(44.667187468730068781 + rm8)/(fac8*r);
            case 9:
                if(r < 1) {
                    double r4 = r2*r2;
                    return 70*(-1339.3818732572790373 + r2)*(-729.14612203399050417 + r2)*
                        (-387.88675827067060694 + r2)*(-182.06453197357443423 + r2)*
                        (247.85684584958031176 - 30.775584314506662372*r2 + r4)*
                        (271.64253396599716354 - 26.122580680557891357*r2 + r4)*
                        (339.79046405460486136 - 14.816743387077397169*r2 + r4)*
                        (567.95216483585341562 + 10.194193917656533542*r2 + r4)/fac9;
                }
                if(r < 3) 
                    return -56*(-34.582311920734983338 + r)*(-24.992661957284416945 + r)*(-17.488028192734710061 + r)*
                        (-1.383150155092837703e-9 + r)*(9.9641140191887944306 + r)*
                        (15.275628751504522331 + r)*(21.325806539922672215 + r)*(28.648730386050029958 + r)*
                        (38.325425823332099091 + r)*(28.562622189135126569 - 10.6719123338422403573*r + r2)*
                        (28.70507115906184515 - 10.5541760339839146887*r + r2)*
                        (29.05913903831001745 - 10.2731288535746689505*r + r2)*
                        (29.910351008219846128 - 9.6992281634742121841*r + r2)*
                        (33.109435948561163527 - 8.5253445394662507222*r + r2)*
                        (15.602088411365139391 + 3.6984662028815353093*r + r2)*
                        (8.74429031620716939 + 4.6043633640202096921*r + r2)*
                        (6.459782693551278648 + 4.944256909578684375*r + r2)/(r*fac9);  
                if(r < 5) 
                    return 28*(-32.102003289043321523 + r)*(-22.256513774013999868 + r)*
                        (-0.019481535609513795422 + r)*(6.7554337530156031483 + r)*
                        (11.394577631317107669 + r)*(16.590631798368255241 + r)*(22.666886952847596752 + r)*
                        (30.066339108752666659 + r)*(39.853988094659807539 + r)*
                        (44.413272765229737826 - 13.3228573032607344464*r + r2)*
                        (44.361461572930625046 - 13.2657621038680272*r + r2)*
                        (44.237151355715103365 - 13.1303141049877793058*r + r2)*
                        (43.97780832824351715 - 12.8542250410841235126*r + r2)*
                        (43.25042168424317674 - 12.2421269588397814766*r + r2)*
                        (43.510305515279639999 - 10.9664711379826025195*r + r2)*
                        (8.706003701431280234 + 1.1855700035587559181*r + r2)*
                        (1.5885655973358157054 + 1.646327906170090721*r + r2)/(r*fac9);
                if(r < 7) {
                    double rm7 = r-7;
                    double rm72 = rm7*rm7;
                    return -8*(-21.671551106850203861 + rm7)*(6.1658351104757744729 + rm7)*
                        (10.795082806936117749 + rm7)*(14.886690765689426014 + rm7)*
                        (19.439779413793422907 + rm7)*(24.671973885987940333 + rm7)*
                        (30.82512647015834664 + rm7)*(38.323916055253040191 + rm7)*
                        (48.234343836705240945 + rm7)*(0.8204755636440703225 - 1.8009337679297521787*rm7 + 
                          rm72)*(0.8875491280110375689 - 1.7846057971090110034*rm7 + rm72)*
                        (1.045466534148576728 - 1.7461307819268706488*rm7 + rm72)*
                        (1.3620925810392051405 - 1.6688498854814113808*rm7 + rm72)*
                        (2.0261760129867730361 - 1.5061616285750895082*rm7 + rm72)*
                        (3.6849852800828341851 - 1.0961867085438506639*rm7 + rm72)*
                        (9.802984897298482732 + 0.4622105022694643882*rm7 + rm72)*
                        (42.534612957672133683 + 12.4694608291474156047*rm7 + rm72)/(r*fac9);
                }
                // (7 < r < 9)
                double rm9 = r - 9;
                double rm9n = rm9*rm9;//(r-9)^2
                rm9n *= rm9n;//(r-9)^4
                rm9n *= rm9n;//(r-9)^8
                rm9n *= rm9n;//(r-9)^16
                return rm9n*(6.825729439553291229 + rm9)*(10.092343822556426859 + rm9)*
                    (13.667628739160225217 + rm9)*(17.697176468096202415 + rm9)*
                    (22.295357761443391034 + rm9)*(27.60630746610677671 + rm9)*
                    (33.854671868419147162 + rm9)*(41.461619200225694147 + rm9)*
                    (51.499165234438845228 + rm9)/(r*fac9);
            case 10:
                if(r < 2) return 126*(-40.730861380535018703 + r)*(-30.704689588645562277 + r)*
                    (-22.992875380046477405 + r)*(-16.360939093305694833 + r)*
                    (12.690092798604652806 + r)*(18.405158447507681799 + r)*(24.824095152789657241 + r)*
                    (32.51599265654231769 + r)*(42.598021872462994313 + r)*
                    (25.300112484166557948 - 10.0334127425641453339*r + r2)*
                    (25.713661830911233875 - 9.8916774102049681353*r + r2)*
                    (26.711889828809740328 - 9.5638283492805680772*r + r2)*
                    (28.888659834912491186 - 8.9292781791424552726*r + r2)*
                    (34.826531125934265278 - 7.6858921844156450972*r + r2)*
                    (25.542386251091209309 + 5.2925509705319046869*r + r2)*
                    (17.55329728692822335 + 6.4205917029794838437*r + r2)*
                    (14.390323965041409498 + 6.9603561206423636862*r + r2)*
                    (13.150208706289443308 + 7.1865945860794790674*r + r2)/fac10;
                if(r < 4) return -84*(-38.554086424835663099 + r)*(-28.493977199105514103 + r)*
                    (-20.519100496194077693 + r)*(-2.8435598074623484647e-6 + r)*
                    (9.3518533371588264352 + r)*(14.376666106461544584 + r)*(19.918636225956956607 + r)*
                    (26.3312575146666767 + r)*(34.076630980093254817 + r)*(44.250454721418781081 + r)*
                    (40.792020913412716658 - 12.76277497587141496*r + r2)*
                    (40.844421107871237868 - 12.6788402943927195074*r + r2)*
                    (40.965626243924349424 - 12.4861948672875134966*r + r2)*
                    (41.210832706836601715 - 12.1196055643902069772*r + r2)*
                    (41.858439634680617615 - 11.4268079110621117711*r + r2)*
                    (45.088363178240802647 - 10.0901947966782361535*r + r2)*
                    (17.28270478888944741 + 2.8693816458246154429*r + r2)*
                    (8.03623013733705846 + 3.8075831701851841313*r + r2)*
                    (4.7572329776597928001 + 4.1491216716114254256*r + r2)/(r*fac10);
                if(r < 6) return 36*(-35.876295340817458023 + r)*(-25.486721357140048453 + r)*
                    (-0.11361367813780578409 + r)*(6.262165656228605773 + r)*
                    (10.733589514338079033 + r)*(15.617956156305644566 + r)*(21.165609509853070804 + r)*
                    (27.633811578402871091 + r)*(35.461618992852383105 + r)*(45.743181407742815902 + r)*
                    (59.115083989577665459 - 15.373313555916510198*r + r2)*
                    (59.055745145993232657 - 15.332428121601498519*r + r2)*
                    (58.919642130832493691 - 15.239153659589212282*r + r2)*
                    (58.65974459256882318 - 15.062822611305076983*r + r2)*
                    (58.16204476708872742 - 14.728077274991556439*r + r2)*
                    (56.77789490034247689 - 13.993336943795713505*r + r2)*
                    (56.5236381230171421 - 12.552528545507928907*r + r2)*
                    (10.404553807463664483 + 0.3594811877983014611*r + r2)*
                    (1.4779285284866008174 + 0.7808770852810373599*r + r2)/(r*fac10);
                if(r < 8) {
                    rm8 = r - 8;
                    rm82 = rm8*rm8;
                    return -9*(-24.174516096409674299 + rm8)*(6.8386234163115958519 + rm8)*
                        (11.380365980370445245 + rm8)*(15.379541701854251628 + rm8)*
                        (19.737652581491112093 + rm8)*(24.638107878947582499 + rm8)*
                        (30.2444382224183403 + rm8)*(36.792247977362783739 + rm8)*
                        (44.715010116238448113 + rm8)*(55.111387208307203333 + rm8)*
                        (0.8333065313057948019 - 1.8173616675857947139*rm8 + rm82)*
                        (0.8864369917459426455 - 1.8056861560167640538*rm8 + rm82)*
                        (1.0071149286455085158 - 1.7791532320571232352*rm8 + rm82)*
                        (1.2333222245296995433 - 1.7293665808953077083*rm8 + rm82)*
                        (1.6563519947115410641 - 1.6360801293141934015*rm8 + rm82)*
                        (2.5208398563960954861 - 1.4447080469488793034*rm8 + rm82)*
                        (4.6566996455212848487 - 0.9676064727655324345*rm8 + rm82)*
                        (12.456219501821072562 + 0.8288654241132394173*rm8 + rm82)*
                        (51.44744926980987549 + 13.688237874578266929*rm8 + rm82)/(r*fac10);
                }
                // 8 < r < 10
                double rm10n = r - 10;// (r-10)
                rm10n = rm10n*rm10n*rm10n;// (r-10)^3
                rm10n = rm10n*rm10n*rm10n;// (r-10)^9
                rm10n = rm10n * rm10n; //(r-10)^18
                return rm10n*(-2.5466026413060609463 + r)*(0.73874700605006930547 + r)*(4.2657055132582012275 + r)*
                    (8.1739240266672909938 + r)*(12.557217578288886452 + r)*(17.519593630574645191 + r)*
                    (23.204770936312821388 + r)*(29.842253660028354113 + r)*(37.865247833846061979 + r)*
                    (48.379142456279730296 + r)/(r*fac10);
            default: throw new IllegalArgumentException(); 
        }
    }

    private static double ArcCoth(double x) {
        return 0.5 * Math.log((x+1)/(x-1));
    }

    public double getChi(double temperature) {
        return 1;
    }

    public void rejectNotify() {
        throw new RuntimeException("nope");
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    protected final double sigma;
    protected final Vector axis0;
    protected final double[] normalWidth;
    protected boolean[] inserted;
    protected long[] numInserts, numTrials;
    protected long totalCount;
    protected final static double fac3 = (9.*(-92 + 27*Math.log(3.)));
    protected final static double fac4 = (4608.*(-179 + 256*ArcCoth(3)));
    protected final static double fac5 = (90.*(-2773712 + 6640625*ArcCoth(4) + 3727*Math.log(3)));
    protected final static double fac6 = -8847360*(-270257 + 905418*ArcCoth(5) + 6245*Math.log(2));
    protected final static double fac7 = -315*(-15825661837824. + 3573372188392.*ArcCoth(4) + 65635383907142.*ArcCoth(6) + 12807215*Math.log(3));
    protected final static double fac8 = -35936796672.*(-1846795758 + 468503*ArcCoth(3) + 1247929879*ArcCoth(5) + 8522825728.*ArcCoth(7));
    protected final static double fac9 = 1.0850962247457076546e23;//-63504*(-10426242710867380224. + 42472111822285718.*ArcCoth(4) + 16104102360845723218.*ArcCoth(6) + 47728241004256868637.*ArcCoth(8) + 1366272181*Math.log(3));
    protected final static double fac10 = 1.2839372233185790761e27;//-72477573120L*(-130183458512569272L + 55413583302L*ArcCoth(3) + 3559737320287686L*ArcCoth(5) + 385043557076036577L*ArcCoth(7) + 505161285400390625L*ArcCoth(9));
    
    public static void main(String[] args) {
//        lensPointTest();
//        System.exit(0);
//        new MCMoveClusterAtomHSRing(new RandomNumberGenerator(), Space3D.getInstance(), 1.0);
//        System.exit(1);
        if (false) {
            double dr = 0.01;
            for (int i=10; i<11; i++) {
                double lasty=0, lastdy=0, lastddy=0, lastd3y = 0;
                for (int jx=1; jx<100*i; jx++) {
                    double r = dr*jx;
                    double y = Math.log(separationProbability(i, r));
                    double dydx = (y-lasty)/dr;
                    double ddydxx = (dydx-lastdy)/dr;
                    double d3y = (ddydxx-lastddy)/dr;
                    double d4y = (d3y-lastd3y)/dr;
                    if (jx>6 && (Double.isInfinite(d3y))) break;
                    if (jx>5) System.out.println(r+" "+Math.abs(y)+" "+Math.abs(dydx)+" "+Math.abs(ddydxx)+" "+Math.abs(d3y)+" "+Math.abs(d4y));
                    lasty = y;
                    lastdy = dydx;
                    lastddy = ddydxx;
                    lastd3y = d3y;
                }
                System.out.println("&");
            }
            System.exit(1);
        }
        double[] b = new double[]{0,1.1,0.9*1.19776,0.674043,0.543464,0.4489,0.381788,0.331708,0.29308,0.262442,0.237574};
        RandomMersenneTwister random = new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray());
        FileWriter fw = null;
//        try {
//            fw = new FileWriter("maxbias.dat");
//        }
//        catch (IOException e) {
//            throw new RuntimeException(e);
//        }
        long maxmaxk = -1;
        for (int i=4; i<=12; i+=2) {
            int j = i/2;
            double normalWidthP = Math.sqrt(0.5/b[j]);
            double normalWidthN = Math.sqrt(0.5/b[i-j]);
            double gaussPrefacP = 1; //2*Math.sqrt(b[j]/Math.PI);
            double gaussPrefacN = 1; //2*Math.sqrt(b[i-j]/Math.PI);
            double v = 0;
            for (int jr=0; jr<100*i; jr++) {
                double rij = jr*0.01;
                double bias = 0;
                double sp = normalWidthP;
                double sp2 = sp*sp;
                double sn = normalWidthN;
                double sn2 = sn*sn;
                double s = Math.sqrt(sp2*sn2/(sp2+sn2));
                double cx = sp2/(sp2+sn2)*rij;
                long nk = 10000000;
                double maxy = -1;
                long maxk = -1;
                for (long k=0; k<=nk; k++) {
                    double x = 0;
                    double jx = x*s;
                    x = k*0.000001;
                    double jy = x*s;
                    // (2*Math.sqrt(b/Math.PI))
                    // (2/(sqrt(b/pi))^3
                    double pGauss = gaussPrefacP*gaussPrefacN * Math.exp(-jy*jy/sp2 - cx*cx/sp2);
                    jx += cx;
                    double rp = Math.sqrt(jx*jx + jy*jy);
                    double rn = Math.sqrt((rij-jx)*(rij-jx) + jy*jy);
                    double pp = separationProbability(j,rp);
                    if (pp == 0) continue;
                    double pn = separationProbability(i-j,rn);
                    if (pn == 0) continue;
                    double p = pp*pn;
                    if (p/pGauss>bias) {
//                        System.out.println(i+" "+rij+" "+(jx-cx)+" "+jy+" "+p/pGauss);
                        bias = p/pGauss;
                        maxy = jy;
                        maxk = k;
                    }
                    else {
//                        break;
                    }
//                    System.out.println(k+" "+jy+" "+p/pGauss);
                }
//                if (true) System.exit(1);
                if (maxy==0) break;
                if (jr==0) v = maxy;
                System.out.println(i+" "+rij+" "+maxy+" "+(v*v - maxy*maxy - 0.25*rij*rij)+" "+bias);
                if (maxk > maxmaxk) maxmaxk = maxk;
                break;
            }
        }
        System.out.println("global max k "+maxmaxk);
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
