package etomica.models.hexane;

import etomica.atom.AtomLeaf;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.normalmode.CoordinateDefinitionMolecule;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;

/**
 * Implementation of CoordinateDefinitionMolecule that handles hexane.
 * CoordinateDefinitionMolecule handles the molecular center of mass. This class
 * adds 3 measures of the molecule orientation. Two that describe the deviation
 * of the molecule's primary orientation and a third that describes rotation
 * about that axis.
 * 
 * @author Andrew Schultz
 * @author Nancy Cribbin
 */
public class CoordinateDefinitionHexane extends CoordinateDefinitionMolecule {

    public CoordinateDefinitionHexane(SpeciesHexane species){
        super(Space3D.getInstance(), 3);

        length = species.getBondLength();
        phi = species.getBondAngle();
        
        axes = new Vector3D[3];
        axes[0] = new Vector3D();
        axes[1] = new Vector3D();
        axes[2] = new Vector3D();
        angles = new Vector3D();
        axis0prime = new Vector3D();
        b = new Vector3D();
        bPrime = new Vector3D();
        midpoint13 = new Vector3D();
        deltaVPrime = new Vector3D();
        c = new Vector3D();
        deltaV = new Vector3D();
        
        vex = new Vector3D();
        temp = new Vector3D();
        axial = new Vector3D();
    }

    public double[] calcU(IAtom molecule) {
        // handle center-of-mass part
        // super.calcU fills in the first 3 elements of |u|
        super.calcU(molecule);

        // Now we play with the molecule we are measuring.

        // Long rotational axis of atom 1
        IVector leafPos1 = ((AtomLeaf) ((IAtomGroup) molecule)
                .getChildList().get(0)).getPosition();
        IVector leafPos2 = ((AtomLeaf) ((IAtomGroup) molecule)
                .getChildList().get(1)).getPosition();
        IVector leafPos3 = ((AtomLeaf) ((IAtomGroup) molecule)
                .getChildList().get(2)).getPosition();
        axis0prime.Ev1Mv2(leafPos3, leafPos1);
        // axis0Prime goes from the 1st atom on the molecule to the 3rd atom on
        // the molecule
        axis0prime.TE(1 / Math.sqrt(axis0prime.squared()));

        // now we do the whole projection thing!! notes 11/3/06 will tell you
        // why this is possible
        // mostly because the axes are unit vectors.
        // Project axis0prime onto axis1 and axis2
        u[3] = axis0prime.dot(axes[1]);
        u[4] = axis0prime.dot(axes[2]);

        // we need to rotate pos2 back to the original frame of reference
        midpoint13.Ev1Pv2(leafPos3, leafPos1);
        midpoint13.TE(0.5);
        // prime = current frame of reference
        // deltaV means the vector coming from the molecule axis and point to
        // the second atom
        deltaVPrime.Ev1Mv2(leafPos2, midpoint13);
        deltaVPrime.TE(1 / Math.sqrt(deltaVPrime.squared()));

        // alpha is the angle between the current molecule orientation and the
        // original
        // molecule orientation
        double cosAlpha = axis0prime.dot(axes[0]);
        if (cosAlpha > 0.999999) {
            // no rotation. if cosAlpha is too close to 1, things blow up, so
            // pretend it's not rotated
            deltaV.E(deltaVPrime);
        } else {
            // bprime is a vector of length 1, normal to axis0prime (the current
            // molecule orientation)
            // and in the plane that contains axis0prime and axes[0] (the
            // original molecule orientation)
            bPrime.Ea1Tv1(1.0 / cosAlpha, axes[0]);
            bPrime.ME(axis0prime);
            bPrime.TE(1 / Math.sqrt(bPrime.squared()));

            // b is a vector of length 1, normal to axis0 (the original molecule
            // orientation) and
            // in the same plane as everything else.
            b.Ea1Tv1(cosAlpha, bPrime);
            // sin^2(alpha) = 1 - cos^2(alpha)
            b.PEa1Tv1(-Math.sqrt(1 - cosAlpha * cosAlpha), axis0prime);

            // calculate the component of deltaVprime in the c direction
            // by subtracting off the component in the bPrime direction
            double bComponent = deltaVPrime.dot(bPrime);

            // get the "c" axis, perpendicular to aprime and bprime (and also a
            // and b)
            c.E(axis0prime);
            c.XE(bPrime);
            double cComponent = deltaVPrime.dot(c);

            // if theta is the angle of rotation of deltaVprime around
            // axis0prime
            // using b as theta=0. We seek to maintain that angle of rotation
            // for
            // deltaV around axis0, using b as theta=0

            deltaV.Ea1Tv1(bComponent, b);
            deltaV.PEa1Tv1(cComponent, c);
        }

        // we want a unit vector. to get the real deltaV, scale it out to
        // deltaVprime's length
        // deltaV.TE(Math.sqrt(deltaVprime.squared()));

        // we don't actually need these
        // newMidpoint13.Ea1Tv1(Math.sqrt(midpoint13.squared()), axes[0]);
        // newPos2.Ev1Pv2(newMidpoint13, deltaV);

        // we just want the angle between the transformed deltaV and axis2
        double dot = deltaV.dot(axes[1]);
        if (dot > 0.9999999) {
            // dot might be >1 due to roundoff, which will make acos blow up
            // also, XE would also return nonsense
            u[5] = 0;
            return u;
        }
        if (dot < -0.99999) {
            // same for -1. This would be really bad as the harmonic potential
            // won't know what to do with this
            // System.out.println("uh-oh, hexane rotated 180");
            u[5] = Math.PI;
            return u;
        }

        u[5] = Math.acos(dot);

        // figure out whether the angle should be negative or not, acos can't
        // tell deltaV x axes[1] will point along the molecule axis or opposite to it
        deltaV.XE(axes[1]);
        // (deltaV x axes[1]) dot axes[0] will be positive if the cross product
        // points
        // in the same direction, but will be negative if the cross product
        // points in
        // the opposite direction
        if (deltaV.dot(axes[0]) < 0) {
            // it's arbitrary which one we flip, so long as we're consistent and
            // do the same thing in setToU
            u[5] = -u[5];
        }
        
      double tol = 0.0000000001;
      for (int i = 0; i < 6 - 3; i++) {
          vex.E(((AtomLeaf) ((IAtomGroup) molecule)
                  .getChildList().get(i)).getPosition());
          vex.ME(((AtomLeaf) ((IAtomGroup) molecule)
                  .getChildList().get(i+1)).getPosition());
          temp.E(((AtomLeaf) ((IAtomGroup) molecule)
                  .getChildList().get(i+3)).getPosition());        
          temp.ME(((AtomLeaf) ((IAtomGroup) molecule)
                  .getChildList().get(i+2)).getPosition());        
          axial.E(((AtomLeaf) ((IAtomGroup) molecule)
                  .getChildList().get(i+2)).getPosition());
          axial.ME(((AtomLeaf) ((IAtomGroup) molecule)
                  .getChildList().get(i+1)).getPosition());
                      
          // Project each vector onto the axial vector, and subtract the
          // axial portion from the result, leaving the radial portion
          axial.normalize();
          axial.TE(length * Math.cos((Math.PI - phi))); // (Pi - phi is
                                                          // 180 - 109.47,
                                                          // or the acute
                                                          // angle)
          vex.ME(axial);
          axial.TE(-1.0); // we do this because we have assumed that
          // the angle between the two vectors (axial
          // & whatever) is obtuse. The angle between
          // axial and temp is acute, so we reverse
          // axial to make the angle acute.
          temp.ME(axial);

          vex.normalize();
          temp.normalize();

          // Calculate the angle between the projected vectors

          /*
           * We have to do some stuff to "fix" answers that are just over
           * 1.0 or under -1.0 because we get NaN, and sometimes we are
           * just really close (like, into the pico and femto range) and
           * we can make that happy.
           */
          double makeGood = vex.dot(temp);
          if (makeGood < -1) {
              if ((makeGood + tol) > -1) {
                  makeGood = -1.0;
              }
          }
          if (makeGood > 1) {
              if ((makeGood - tol) < 1) {
                  makeGood = 1.0;
              }
          }
          
          u[6+i] = Math.acos(makeGood) - angles.x(i);
          
      }//end of the for loop that calculates the angles and stores the differences.

        
        
        
        
        
        return u;
    }

    public void initNominalU(IAtom molecule) {
        // handle center-of-mass part
        super.initNominalU(molecule);
        // assume they're all oriented the same way.
        
        // Set up all the axes based on the molecule atom0, the reference
        // molecule
        // Long rotational axis of atom 0
        IVector leafPos1 = ((AtomLeaf) ((IAtomGroup) molecule)
                .getChildList().get(0)).getPosition();
        IVector leafPos2 = ((AtomLeaf) ((IAtomGroup) molecule)
                .getChildList().get(1)).getPosition();
        IVector leafPos3 = ((AtomLeaf) ((IAtomGroup) molecule)
                .getChildList().get(2)).getPosition();
        // axes[0] should point from the 0th atom on the molecule to the 2nd
        // atom on the molecule
        axes[0].Ev1Mv2(leafPos3, leafPos1);
        // Now we take the midpoint between the 0th atom and the 2nd atom.
        axes[1].Ev1Pv2(leafPos1, leafPos3);
        axes[1].TE(-0.5);
        // Then we subtract the midpoint from the location of the 1st atom
        // on the molecule to get our final vector.
        // AKA isosceleshappyvector
        axes[1].PE(leafPos2);

        // Normalize our axes
        axes[0].TE(1.0 / Math.sqrt(axes[0].squared()));
        axes[1].TE(1.0 / Math.sqrt(axes[1].squared()));

        // Last axis is simply normal to the other two
        axes[2].E(axes[0]);
        axes[2].XE(axes[1]);
        
//        double tol = 0.0000000001;
//        for (int i = 0; i < 6 - 3; i++) {
//            vex.E(((AtomLeaf) ((IAtomGroup) molecule)
//                    .getChildList().get(i)).getPosition());
//            vex.ME(((AtomLeaf) ((IAtomGroup) molecule)
//                    .getChildList().get(i+1)).getPosition());
//            temp.E(((AtomLeaf) ((IAtomGroup) molecule)
//                    .getChildList().get(i+3)).getPosition());        
//            temp.ME(((AtomLeaf) ((IAtomGroup) molecule)
//                    .getChildList().get(i+2)).getPosition());        
//            axial.E(((AtomLeaf) ((IAtomGroup) molecule)
//                    .getChildList().get(i+2)).getPosition());
//            axial.ME(((AtomLeaf) ((IAtomGroup) molecule)
//                    .getChildList().get(i+1)).getPosition());
//                        
//            // Project each vector onto the axial vector, and subtract the
//            // axial portion from the result, leaving the radial portion
//            axial.normalize();
//            axial.TE(length * Math.cos((Math.PI - phi))); // (Pi - phi is
//                                                            // 180 - 109.47,
//                                                            // or the acute
//                                                            // angle)
//            vex.ME(axial);
//            axial.TE(-1.0); // we do this because we have assumed that
//            // the angle between the two vectors (axial
//            // & whatever) is obtuse. The angle between
//            // axial and temp is acute, so we reverse
//            // axial to make the angle acute.
//            temp.ME(axial);
//
//            vex.normalize();
//            temp.normalize();
//
//            // Calculate the angle between the projected vectors
//
//            /*
//             * We have to do some stuff to "fix" answers that are just over
//             * 1.0 or under -1.0 because we get NaN, and sometimes we are
//             * just really close (like, into the pico and femto range) and
//             * we can make that happy.
//             */
//            double makeGood = vex.dot(temp);
//            if (makeGood < -1) {
//                if ((makeGood + tol) > -1) {
//                    makeGood = -1.0;
//                }
//            }
//            if (makeGood > 1) {
//                if ((makeGood - tol) < 1) {
//                    makeGood = 1.0;
//                }
//            }
//            angles.setX(i, Math.acos(makeGood));
//        }//end of the for loop that calculate the initial angles and stores them.
        
        //We know this from the initial model; the intial torsional angles are PI
        for(int i = 0; i < 3; i++){
            angles.setX(i, Math.PI);
        }
    }

    public void setToU(IAtom atom, double[] u) {
        super.setToU(atom, u);
        throw new RuntimeException("Don't yet know how to set orientation");
    }

    private static final long serialVersionUID = 1L;

    private final Vector3D axis0prime, bPrime, midpoint13, deltaVPrime;

    private final Vector3D c, deltaV, b;

    private Vector3D[] axes;
    
    private Vector3D angles;
    
    private Vector3D vex, temp, axial;
    
    private double length, phi;

    public void setLength(double length) {
        this.length = length;
    }

    public void setPhi(double phi) {
        this.phi = phi;
    }
}