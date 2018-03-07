/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.exception.MethodNotImplementedException;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePositionCOM;
import etomica.molecule.OrientationCalc;
import etomica.molecule.OrientationCalcQuaternion;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.RotationTensor3D;

public class OrientationCalcWater3P extends ConformationWater3P implements
        OrientationCalc, OrientationCalcQuaternion, java.io.Serializable {

    public OrientationCalcWater3P(Space space) {
        super(space);
        xWork = space.makeVector();
        yWork = space.makeVector();
        zWork = space.makeVector();
        com0 = space.makeVector();
        rotationTensor = (RotationTensor3D)space.makeRotationTensor();
        atomPositionCOM = new MoleculePositionCOM(space);
        initialized = false;
        previousTensor = space.makeRotationTensor();
        workTensor = space.makeRotationTensor();
        double hMass = 1.0079;
        double oMass = 15.9994;
        com0.E(new double[] {hMass*bondLengthOH, 0, 0.0});
        com0.PEa1Tv1(hMass, space.makeVector(new double[] {bondLengthOH*Math.cos(angleHOH), bondLengthOH*Math.sin(angleHOH), 0.0}));
        com0.TE(1.0/(hMass*2+oMass));
    }

    public void calcOrientation(IMolecule molecule, double[] quat) {
        // depend on ordering as H1, H2, O.  we could sniff this out if needed
        IAtomList children = molecule.getChildList();
        IAtom H1 = children.get(0);
        IAtom H2 = children.get(1);
        IAtom O = children.get(2);

        // xyzWork define the orientation of the molecule
        // xWork is the axis defined by the H1->H2 vector
        // yWork is perpendicular to xWork, from O to the point midway between H1 and H2
        // zWork is perpendicular to xWork and yWork
        
        xWork.Ev1Mv2(H2.getPosition(), H1.getPosition());
        xWork.normalize();

        yWork.Ev1Pv2(H1.getPosition(), H2.getPosition());
        yWork.TE(0.5);
        yWork.ME(O.getPosition());
        yWork.normalize();

        zWork.E(xWork);
        zWork.XE(yWork);
        
        double xx = xWork.getX(0);
//        double yx = yWork.x(0);
        double zx = zWork.getX(0);
        double xy = xWork.getX(1);
//        double yy = yWork.x(1);
        double zy = zWork.getX(1);
        // azy = |zy|, which appears in expressions for cos(psi) 
        double azy = Math.abs(zy);
//        double yz = yWork.x(2);
        double xz = xWork.getX(2);
        double zz = zWork.getX(2);
        // fix roundoff errors
        if (zz > 1) zz = 1;
        else if (zz < -1) zz = -1;
        // alpha is just a term that appears in the trig functions of psi and quaternion expressions
        double alpha = Math.sqrt(zx*zx+zy*zy);
        
//        System.out.println(xx+" "+xy+" "+xz+" "+yx+" "+yy+" "+yz+" "+zx+" "+zy+" "+zz+" "+alpha);
//        System.out.println(xx+" "+xy+" "+xz+" yx yy yz "+zx+" "+zy+" "+zz+" "+alpha);
        
        
        // I hate Euler angles
        // theta = acos(zz)
        // psi = atan(-xz/yz)
        // phi = acos(xx cos(psi) + yx sin(psi))
        //     = acos(xx*yz/alpha - yx*xz/alpha)
        // psi the the last rotation applied.
        // psi is a rotation about the original z axis.  positive psi means
        // x is rotated into the direction of y.  -pi/2 < psi < +pi/2
        
        // theta is the second rotation applied.
        // theta is a rotation about the original x axis.  positive theta
        // means y is rotated into the direction of z.  -pi < psi < +pi
        
        // phi is the first rotation applied.
        // phi is a rotation about the original z axis.  positive phi means
        // x is rotated into the direction of y.  -pi < phi < +pi
        
        boolean psiNegative = zx*zy > 0;
        boolean thetaNegative = zy > 0;
        boolean phiNegative = xz < 0 != thetaNegative;

        double szx = Math.abs(zx);
        if (psiNegative) {
            szx = -szx;
            //we'll use sin(psi) = sxz / alpha
            // sin(psi) < 0 if psi < 0, so we must impose that on our expression
            // regardless of the sign of xz
        }
        double costheta = zz;
        double sintheta = Math.sqrt(1-zz*zz);
        if (thetaNegative) {
            sintheta = -sintheta;
        }
        double sinphi = xz/sintheta;// assumes -pi/2 < phi < pi
        // fix roundoff errors
        if (sinphi > 1) sinphi = 1;
        else if (sinphi < -1) sinphi = -1;
            
        double cosphi = Math.sqrt(1.0-sinphi*sinphi);

        double cospsi = azy / alpha;
        double sinpsi = szx / alpha;
        
//        double theta = Math.acos(zz);
//        if (thetaNegative) theta = -theta;
//        double psi = -Math.atan(zx/zy);
//        double phi = Math.asin(xz/Math.sin(theta));
//        System.out.println("phi, theta, psi "+phi+" "+theta+" "+psi);
//        System.out.println("costheta sintheta "+costheta+" "+sintheta+" "+Math.cos(theta)+" "+Math.sin(theta));
//        System.out.println("cosphi sinphi "+cosphi+" "+sinphi+" "+Math.cos(phi)+" "+Math.sin(phi));
//        System.out.println("cospsi sinpsi "+cospsi+" "+sinpsi+" "+Math.cos(psi)+" "+Math.sin(psi));
//        phi = Math.acos(xx*ayz/alpha + yx*sxz/alpha);
//        System.out.println("phi2 "+phi);
        
        // q0 = cos(0.5theta)cos(0.5(phi+psi))
        // q1 = sin(0.5theta)cos(0.5(phi-psi))
        // q2 = sin(0.5theta)sin(0.5(phi-psi))
        // q3 = cos(0.5theta)sin(0.5(phi+psi))

        // cos(0.5(A+B)) = cos(0.5A+0.5B) = cos(0.5A)cos(0.5B)-sin(0.5A)sin(0.5B)
        // sin(0.5(A+B)) = sin(0.5A+0.5B) = sin(0.5A)cos(0.5B)+cos(0.5A)sin(0.5B)
        // cos(0.5A) = sqrt((1+cosA)/2)
        // sin(0.5A) = Z sqrt((1-cosA)/2) (Z=1 if 0<0.5A<pi, Z=-1 if -pi<0.5A<0)

        // these are actually half-angle (cos(0.5phi)*cos(0.5*psi))
        double cosPhicosPsi;
        double sinPhisinPsi;
        double cosPhisinPsi;
        double sinPhicosPsi;
        if (Math.abs(zz) > 0.99999999) {
            // theta is very small, so our calculation of phi and psi would be grossly imprecise 
//            theta = 0;
            sintheta = 0;
            costheta = 1;
//            psi = 0;
            // -pi < phi < +pi.  oops
//            phi = Math.acos(xx);
//            if (xy < 0) phi = -phi;
//            System.out.println("new phi, theta, psi "+phi+" "+theta+" "+psi);

            // we would take psi to be atan(-xz/yz), which fails
            // take psi to be 0 and let phi handle rotations
            sinPhisinPsi = 0;
            cosPhisinPsi = 0;
            
            cosPhicosPsi = Math.sqrt((1.0+xx)/2.0);
            sinPhicosPsi = Math.sqrt((1.0-xx)/2.0);

            psiNegative = false;
            phiNegative = xy < 0;
            thetaNegative = false;
        }
        else {
            if (xx*cospsi + xy*sinpsi < 0) {
                // this doesn't work well if x is almost pointing in the z direction
                // but that would correspond to phi = pi/2, so getting this part right
                // doesn't really matter (cosphi ~= 0)
//                if (phi > 0) {
//                    phi = Math.PI - phi;
//                }
//                else {
//                    phi = -Math.PI - phi;
//                }
//                System.out.println("new phi "+phi);
                cosphi = -cosphi;
            }
//            System.out.println("new phi "+phi);
//            System.out.println("new cosphi "+cosphi);
            cosPhicosPsi = 0.5*Math.sqrt((1.0+cosphi)*(1.0+cospsi));
            sinPhisinPsi = 0.5*Math.sqrt((1.0-cosphi)*(1.0-cospsi));
            cosPhisinPsi = 0.5*Math.sqrt((1.0+cosphi)*(1.0-cospsi));
            sinPhicosPsi = 0.5*Math.sqrt((1.0-cosphi)*(1.0+cospsi));
        }
//        System.out.println("cosphi cospsi "+cosPhicosPsi+" "+(Math.cos(0.5*phi)*Math.cos(0.5*psi)));
//        System.out.println("sinphi sinpsi "+sinPhisinPsi+" "+(-Math.sin(0.5*phi)*Math.sin(0.5*psi)));
//        System.out.println("sinphi "+Math.sqrt((yx*xz-xx*yz+alpha)/(2.0*alpha))+" "+(Math.sin(0.5*phi)));
//        System.out.println("sinpsi "+Math.sqrt((alpha-yz)/(2.0*alpha))+" "+(Math.sin(0.5*psi)));
//        System.out.println("sinpsi "+xz/alpha+" "+(Math.sin(psi)));
//        System.out.println(cosPhicosPsi+" "+sinPhisinPsi+" "+cosPhisinPsi+" "+sinPhicosPsi);
        if (psiNegative) {
            sinPhisinPsi = -sinPhisinPsi;
            cosPhisinPsi = -cosPhisinPsi;
        }
        if (phiNegative) {
            sinPhisinPsi = -sinPhisinPsi;
            sinPhicosPsi = -sinPhicosPsi;
        }
//        System.out.println("cosphi cospsi "+cosPhicosPsi+" "+(Math.cos(0.5*phi)*Math.cos(0.5*psi)));
//        System.out.println("sinphi sinpsi "+sinPhisinPsi+" "+(-Math.sin(0.5*phi)*Math.sin(0.5*psi)));
//        System.out.println("sinphi "+Math.sqrt((alpha-yx*xz-xx*yz)/(2.0*alpha))+" "+(Math.sin(0.5*phi)));
//        System.out.println("sinpsi "+Math.sqrt((alpha-ayz)/(2.0*alpha))+" "+(Math.sin(0.5*psi)));
//        System.out.println("sinpsi "+xz/alpha+" "+(Math.sin(psi)));
//        System.out.println("cosphi sinpsi "+cosPhisinPsi+" "+(Math.cos(0.5*phi)*Math.sin(0.5*psi)));
//        System.out.println("sinphi cospsi "+cosPhisinPsi+" "+(Math.sin(0.5*phi)*Math.cos(0.5*psi)));
//        System.out.println("cosphi "+Math.sqrt((xx*yz+yx*xz+alpha)/(2.0*alpha))+" "+(Math.cos(0.5*phi)));
//        System.out.println("cospsi "+Math.sqrt((ayz+alpha)/(2.0*alpha))+" "+(Math.cos(0.5*phi)));
//        System.out.println("halves "+cosPhicosPsi+" "+sinPhisinPsi+" "+cosPhisinPsi+" "+sinPhicosPsi);
//        System.out.println("halves "+Math.cos(0.5*phi)*Math.cos(0.5*psi)+" "+Math.sin(0.5*phi)*Math.sin(0.5*psi)
//                +" "+Math.cos(0.5*phi)*Math.sin(0.5*psi)+" "+Math.sin(0.5*phi)*Math.cos(0.5*psi));
        double cosHTheta = Math.sqrt((1.0+costheta)/2.0);
        double sinHTheta = Math.sqrt((1.0-costheta)/2.0);
        if (thetaNegative) {
            sinHTheta = -sinHTheta;
        }
//        System.out.println("half thetas "+sinHTheta+" "+cosHTheta);
//        System.out.println("half thetas "+Math.sin(0.5*theta)+" "+Math.cos(0.5*theta));
        quat[0] = cosHTheta*(cosPhicosPsi - sinPhisinPsi);
//        System.out.println("q0 "+quat[0]+" "+Math.cos(0.5*theta)*Math.cos(0.5*(phi+psi)));
        quat[1] = sinHTheta*(cosPhicosPsi + sinPhisinPsi);
//        System.out.println("q1 "+quat[1]+" "+Math.sin(0.5*theta)*Math.cos(0.5*(phi-phi)));
        quat[2] = sinHTheta*(-sinPhicosPsi + cosPhisinPsi);
//        System.out.println("q2 "+quat[2]+" "+Math.sin(0.5*theta)*Math.sin(0.5*(psi-phi)));
        quat[3] = cosHTheta*(sinPhicosPsi + cosPhisinPsi); 
//        System.out.println("q3 "+quat[3]+" "+Math.cos(0.5*theta)*Math.sin(0.5*(phi+psi)));
//        System.out.println((quat[0]*quat[0]+quat[1]*quat[1]+quat[2]*quat[2]+quat[3]*quat[3]));
    }

    public void setOrientation(IMolecule molecule, double[] quaternions) {
        throw new MethodNotImplementedException("oops");
    }

    public void initializePositions(IAtomList childList) {
        super.initializePositions(childList);
        if (!initialized) {
            com0.E(atomPositionCOM.position(childList.get(0).getParentGroup()));
            initialized = true;
        }
    }
    
    private static final long serialVersionUID = 1L;

    protected final Vector xWork, yWork, zWork;
    protected final Vector com0;
    protected final RotationTensor3D rotationTensor;
    protected final MoleculePositionCOM atomPositionCOM;
    protected boolean initialized;
    protected final RotationTensor previousTensor, workTensor;
    
    protected static void doTransform(IMolecule molecule, Vector r0, RotationTensor rotationTensor) {
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }

    public void calcOrientation(IMolecule molecule,
            IOrientationFull3D orientation) {
        // depend on ordering as H1, H2, O.  we could sniff this out if needed
        IAtomList children = molecule.getChildList();
        IAtom H1 = children.get(0);
        IAtom H2 = children.get(1);
        IAtom O = children.get(2);

        // xyzWork define the orientation of the molecule
        // xWork is the axis defined by the H1->H2 vector
        // yWork is perpendicular to xWork, from O to the point midway between H1 and H2
        
        xWork.Ev1Mv2(H2.getPosition(), H1.getPosition());
//        xWork.normalize();

        yWork.Ev1Pv2(H1.getPosition(), H2.getPosition());
        yWork.TE(0.5);
        yWork.ME(O.getPosition());
//        yWork.normalize();

        orientation.setDirections(xWork, yWork);
    }

    public void setOrientation(IMolecule molecule,
            IOrientationFull3D orientation) {
        xWork.E(atomPositionCOM.position(molecule));
        IAtomList childList = molecule.getChildList();
        initializePositions(childList);
        rotationTensor.setOrientation(orientation);
        rotationTensor.invert();
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(com0);
            rotationTensor.transform(r);
            r.PE(xWork);
        }
    }
}
