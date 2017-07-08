/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Hard core with an attractive tail that goes to zero linearly with r.
 * 
 * @author Jhumpa Adhikari
 */

public class P2TriangleWell extends Potential2 {

    public P2TriangleWell(Space space) {
        this(space, 1.0, 2.0, 1.0);
    }
  
    public P2TriangleWell(Space space, double coreDiameter, double lambda, double epsilon) {
        super(space);
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        force = space.makeVector();
        dr = space.makeVector();
    }

    public double energy(IAtomList pair) {
        IAtom atom0 = pair.getAtom(0);
        IAtom atom1 = pair.getAtom(1);
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        boundary.nearestImage(dr);

        double r2 = dr.squared();
       
        if(r2 < coreDiameterSquared)
            return Double.POSITIVE_INFINITY;

        if(r2 > wellDiameterSquared)
            return 0.0;
            
        double r1 = Math.sqrt(r2);
            
        return (epsilon/(lambda - 1.0))*((r1/coreDiameter)- lambda);
    }
 

    // what could call this?
    public Vector force(IAtomList pair){
        
        IAtom atom0 = pair.getAtom(0);
        IAtom atom1 = pair.getAtom(1);
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        if(r2 > wellDiameterSquared){
            force.E(0.0);
        }
        else {
            force.E(dr);
            force.TE(constant/Math.sqrt(r2));//lambda > 1.0
        }
        return force;
    }
     
    public double getCoreDiameter() {return coreDiameter;}
    public final void setCoreDiameter(double c) {
        coreDiameter = c;
        coreDiameterSquared = c*c;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
        constant = epsilon/(coreDiameter*(1.0 - lambda));
    }
    public final Dimension getCoreDiameterDimension() {
        return Length.DIMENSION;
    }
    
    public double getRange() {
        return wellDiameter;
    }

    public double getLambda() {return lambda;}
    public final void setLambda(double lam) {
        lambda = lam;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
        constant = epsilon/(coreDiameter*(1.0 - lambda));
    }
    public final Dimension getLambdaDimension() {
        return Null.DIMENSION;
    }

    public double getEpsilon() {return epsilon;}
    public final void setEpsilon(double eps) {
        epsilon = eps;
        constant = epsilon/(coreDiameter*(1.0 - lambda));
    }
    public final Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    private static final long serialVersionUID = 1L;
    private double coreDiameter, coreDiameterSquared;
    private double wellDiameter, wellDiameterSquared;
    private double lambda; //wellDiameter = coreDiameter * lambda ;lambda is well width
    private double epsilon;
    private double constant;
    private final Vector force;
    private final Vector dr;
    private Boundary boundary;
}

