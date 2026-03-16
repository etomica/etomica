package etomica.starpolymer;

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
    protected double sigma = 1.0, sigma0 = 1.15;
    protected Vector[] graftVectors;

    public ConformationStarPolymerGraft(Space space, int f, int l) {
        this.space = space;
        this.f = f;
        this.l = l;
        makeGraftVectors();
    }

    public void setSigma(double sigma) {
        this.sigma = sigma;
    }

    public void setSigma0(double sigma0) {
        this.sigma0 = sigma0;
    }
    private  void  makeGraftVectors( ) {
        graftVectors = space.makeVectorArray(f);
        if (f == 1) {
            graftVectors[0].setX(0,1);
        }
        else if (f == 2) {
            graftVectors[0].setX(0,1); graftVectors[1].setX(0,-1);
        }
        else if (f == 4) { // tetrahedron
            double x = 1.0/Math.sqrt(3);
            graftVectors[0].E(new double[]{ x, x, x});
            graftVectors[1].E(new double[]{-x,-x, x});
            graftVectors[2].E(new double[]{-x, x,-x});
            graftVectors[3].E(new double[]{ x,-x,-x});
        }
        else if (f == 6) { // octahedron
            graftVectors[0].setX(0,1);  graftVectors[1].setX(0,-1);
            graftVectors[2].setX(1,1);  graftVectors[3].setX(1,-1);
            graftVectors[4].setX(2,1);  graftVectors[5].setX(2,-1);
        }
        else if (f == 8) { // cube corners
            Vector v = space.makeVector(); v.E(1);
            for (int i=0; i<8; i++) {
                graftVectors[i].E(1.0/Math.sqrt(3)); graftVectors[i].TE(v);
                for (int j=2; j>-1; j--) { if (v.getX(j)==1){ v.setX(j,-1); break; } v.setX(j,1); }
            }
        }
        else if (f == 12) { // 12 “FCC-like”
            double phi=(1+Math.sqrt(5))/2, x1=1/Math.sqrt(phi*phi+1), x2=phi*x1;
            graftVectors[0].E(new double[]{  0, x1, x2});
            graftVectors[1].E(new double[]{  0, x1,-x2});
            graftVectors[2].E(new double[]{  0,-x1, x2});
            graftVectors[3].E(new double[]{  0,-x1,-x2});
            graftVectors[4].E(new double[]{ x1, x2,  0});
            graftVectors[5].E(new double[]{ x1,-x2,  0});
            graftVectors[6].E(new double[]{-x1, x2,  0});
            graftVectors[7].E(new double[]{-x1,-x2,  0});
            graftVectors[8].E(new double[]{ x2,  0, x1});
            graftVectors[9].E(new double[]{ x2,  0,-x1});
            graftVectors[10].E(new double[]{-x2,  0, x1});
            graftVectors[11].E(new double[]{-x2,  0,-x1});
        }
        else if (f == 20) { // 20-point set (dodeca/icosa related)
            double phi=(1+Math.sqrt(5))/2, x1=1/Math.sqrt(3), x2=x1*phi, x3=x1/phi;
            graftVectors[0].E(new double[]{ x1, x1, x1}); graftVectors[1].E(new double[]{ x1, x1,-x1});
            graftVectors[2].E(new double[]{ x1,-x1, x1}); graftVectors[3].E(new double[]{ x1,-x1,-x1});
            graftVectors[4].E(new double[]{-x1, x1, x1}); graftVectors[5].E(new double[]{-x1, x1,-x1});
            graftVectors[6].E(new double[]{-x1,-x1, x1}); graftVectors[7].E(new double[]{-x1,-x1,-x1});
            graftVectors[8].E(new double[]{  0, x3, x2}); graftVectors[9].E(new double[]{  0, x3,-x2});
            graftVectors[10].E(new double[]{ 0,-x3, x2}); graftVectors[11].E(new double[]{ 0,-x3,-x2});
            graftVectors[12].E(new double[]{ x3, x2,  0}); graftVectors[13].E(new double[]{ x3,-x2,  0});
            graftVectors[14].E(new double[]{-x3, x2,  0}); graftVectors[15].E(new double[]{-x3,-x2,  0});
            graftVectors[16].E(new double[]{ x2,  0, x3}); graftVectors[17].E(new double[]{ x2,  0,-x3});
            graftVectors[18].E(new double[]{-x2,  0, x3}); graftVectors[19].E(new double[]{-x2,  0,-x3});
        }
        else {
            throw new RuntimeException("Unsupported f="+f+" (use 1,2,4,6,8,12,20)");
        }

        // normalize
        for (int i=0; i<f; i++) {
            double s2 = graftVectors[i].squared();
            graftVectors[i].TE(1.0/Math.sqrt(s2));
        }

    }

    @Override
    public void initializePositions(IAtomList atomList) {
        IAtom core = atomList.get(0);
        core.getPosition().E(0);
        Vector corePos = core.getPosition();

        int iBead = 1;
        Vector dr = space.makeVector();
        Vector lastPos = space.makeVector();

        for (int j = 0; j < f; j++) {
            lastPos.E(corePos);
            for (int k = 0; k < l; k++) {
                IAtom atom = atomList.get(iBead);
                dr.E(graftVectors[j]);

                Vector p = atom.getPosition();
                p.E(lastPos);
                double r01 = 0.5*(sigma0 + sigma);// --new
                p.PEa1Tv1(k == 0 ? r01 : sigma, dr);
                iBead++;
                lastPos.E(p);
            }
        }

    }
}
