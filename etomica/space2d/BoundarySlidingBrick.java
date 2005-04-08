package etomica.space2d;

import etomica.Phase;
import etomica.integrator.IntegratorMD;
import etomica.space.BoundaryRectangular;
import etomica.space.Coordinate;
import etomica.space.Vector;

/**
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Style - Code Templates
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jan 27, 2005 by kofke
 */
public class BoundarySlidingBrick extends BoundaryRectangular {

    private double gamma = 0.0;
    private double delvx;

    private IntegratorMD.ChronoMeter timer;

    public BoundarySlidingBrick() {
        super();
    }

    public BoundarySlidingBrick(Phase p) {
        super(p);
    }

    public void setShearRate(double g) {
        gamma = g;
        computeDelvx();
    }

    public double getShearRate() {
        return gamma;
    }

    private void computeDelvx() {
        delvx = gamma * dimensions.y;
    }

    public void setTimer(IntegratorMD.ChronoMeter t) {
        timer = t;
    }

    public void nearestImage(Vector dr) {
        double delrx = delvx * timer.currentValue();
        double cory = ((dr.y > 0.0) ? Math.floor(dr.y / dimensions.y + 0.5)
                : Math.ceil(dr.y / dimensions.y - 0.5));
        dr.x -= cory * delrx;
        dr.x -= dimensions.x
                * ((dr.x > 0.0) ? Math.floor(dr.x / dimensions.x + 0.5) : Math
                        .ceil(dr.x / dimensions.x - 0.5));
        dr.y -= dimensions.y * cory;
    }

    public void centralImage(Vector r) {
        double delrx = delvx * timer.currentValue();
        double cory = ((r.y >= 0.0) ? Math.floor(r.y / dimensions.y) : Math
                .ceil(r.y / dimensions.y - 1.0));
        //        if(cory != 0.0) System.out.println(delrx*cory);
        r.x -= cory * delrx;
        r.x -= dimensions.x
                * ((r.x >= 0.0) ? Math.floor(r.x / dimensions.x) : Math
                        .ceil(r.x / dimensions.x - 1.0));
        r.y -= dimensions.y * cory;
    }

    public void centralImage(Coordinate c) {
        Vector r = c.r;
        double cory = ((r.y > 0.0) ? Math.floor(r.y / dimensions.y) : Math
                .ceil(r.y / dimensions.y - 1.0));
        double corx = ((r.x > 0.0) ? Math.floor(r.x / dimensions.x) : Math
                .ceil(r.x / dimensions.x - 1.0));
        if (corx == 0.0 && cory == 0.0)
            return;
        double delrx = delvx * timer.currentValue();
        Vector p = c.p;
        r.x -= cory * delrx;
        r.x -= dimensions.x * corx;
        r.y -= dimensions.y * cory;
        p.x -= cory * delvx;
    }

    public double[][] imageOrigins(int nShells) {
        int nImages = (2 * nShells + 1) * (2 * nShells + 1) - 1;
        double[][] origins = new double[nImages][D];
        int k = 0;
        for (int i = -nShells; i <= nShells; i++) {
            for (int j = -nShells; j <= nShells; j++) {
                if (i == 0 && j == 0) {
                    continue;
                }
                origins[k][0] = i * dimensions.x + j * delvx
                        * timer.currentValue();
                origins[k][1] = j * dimensions.y;
                k++;
            }
        }
        return origins;
    }//end of imageOrigins
}//end of BoundarySlidingBrick
