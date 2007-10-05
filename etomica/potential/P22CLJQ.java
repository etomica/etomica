package etomica.potential;

import etomica.atom.AtomSet;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.space.IVector;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;

/**
 * Two-centered Lennard Jones molecule with a quadrupole.
 *
 * @author Jayant K. Singh
 */
public class P22CLJQ extends Potential2 {

    public P22CLJQ(Space space) {
        this(space, 1, 1, 1);
    }

    public P22CLJQ(Space space, double sigma, double epsilon,  double moment) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setQuadrupolarMomentSquare(moment);
        com1 = space.makeVector();
        com2 = space.makeVector();
        v12 = space.makeVector();
        tvect = space.makeVector();
        v1 = space.makeVector();
        v2 = space.makeVector();
        az1 = space.makeVector();
        az2 = space.makeVector();
        dr = space.makeVector();
    }

    public void setHardCoreDiamterSq(double val){
        hsdiasq=val;
    }

    public void setBox(Box box) {
        nearestImageTransformer = box.getBoundary();
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(AtomSet pair){
        double ener=0.0;

        AtomSet mol1 = ((IAtomGroup)pair.getAtom(0)).getChildList(); 
        AtomSet mol2 = ((IAtomGroup)pair.getAtom(1)).getChildList(); 
        IAtomPositioned bead11 = (IAtomPositioned)mol1.getAtom(0);
        IAtomPositioned bead12 = (IAtomPositioned)mol1.getAtom(1);

        IAtomPositioned bead21 = (IAtomPositioned)mol2.getAtom(0);
        IAtomPositioned bead22 = (IAtomPositioned)mol2.getAtom(1);

        // LJ contributation

        dr.Ev1Mv2(bead11.getPosition(), bead21.getPosition());
        nearestImageTransformer.nearestImage(dr);
        ener=ener+calculateEnergy(dr.squared());

        dr.Ev1Mv2(bead11.getPosition(), bead22.getPosition());
        nearestImageTransformer.nearestImage(dr);
        ener=ener+calculateEnergy(dr.squared());

        dr.Ev1Mv2(bead12.getPosition(), bead21.getPosition());
        nearestImageTransformer.nearestImage(dr);
        ener=ener+calculateEnergy(dr.squared());

        dr.Ev1Mv2(bead12.getPosition(), bead22.getPosition());
        nearestImageTransformer.nearestImage(dr);
        ener=ener+calculateEnergy(dr.squared());

        if(Q2!=0.0 && !Double.isInfinite(ener)){

            com1.Ev1Pv2(bead11.getPosition(), bead12.getPosition());
            com1.TE(0.5);

            com2.Ev1Pv2(bead21.getPosition(), bead22.getPosition());
            com2.TE(0.5);
            v12.Ev1Mv2(com1, com2);
            nearestImageTransformer.nearestImage(v12);
            double r2=v12.squared();

            v12.normalize();

            initializeRMAboutRandomVector(v12);
            convertToBodyFrame(v12);

            // axis one
            dr.Ev1Mv2(bead11.getPosition(), bead12.getPosition());
            nearestImageTransformer.nearestImage(dr);

            v1.E(dr);

            convertToBodyFrame(v1);
            v1.normalize();

            //axis two
            dr.Ev1Mv2(bead21.getPosition(), bead22.getPosition());
            nearestImageTransformer.nearestImage(dr);

            v2.E(dr);

            convertToBodyFrame(v2);
            v2.normalize();

            // theta's are the angle between the axis and center -center of the molecule
            double theta1=Math.acos(v1.dot(v12));
            double theta2=Math.acos(v2.dot(v12));

            double c1=Math.cos(theta1);
            double c2=Math.cos(theta2);
            double s1=Math.sin(theta1);
            double s2=Math.sin(theta2);
            az1.Ea1Tv1(s1,v1);
            az2.Ea1Tv1(s2,v2);

            az1.setX(2,0.0);
            az2.setX(2,0.0);
            az1.normalize();
            az2.normalize();

            double c12 = az1.dot(az2);

            double uqq=(3.0*Q2)/(4.0*r2*r2*Math.sqrt(r2));
            double temp1=s1*s2*c12-4*c1*c2;
            ener += uqq*(1.0-5.0*(c1*c1+c2*c2+3*c1*c1*c2*c2) +2*(temp1*temp1));
        }
        return ener;
    }

    protected double calculateEnergy(double r2){
        if(r2<hsdiasq*sigma2) {
            return Double.POSITIVE_INFINITY;
        }
        double s2 = sigma2/(r2);
        double s6 = s2*s2*s2;
        return epsilon4*s6*(s6 - 1.0);
    }

    protected void initializeRMAboutRandomVector(IVector vector1){
        //  vector1 is already normalized
        // initialize WorldUp
        worldUp[0] = 0;
        worldUp[1] = 1.0;
        worldUp[2] = 0;

        // first, calculate and normalize the view vector
        viewOut[0] = vector1.x(0);
        viewOut[1] = vector1.x(1);
        viewOut[2] = vector1.x(2);

        double viewMagnitude=Math.sqrt(vector1.squared());

        // invalid points (not far enough apart)
        if (viewMagnitude < .000001) {
            throw new RuntimeException("invalid points, not far enough apart");
        }
    
        // normalize. This is the unit vector in the "new Z" direction

        // Now the hard part: The ViewUp or "new Y" vector
        double upProjection = viewOut[1];

        viewUp[0] =  - upProjection*viewOut[0];
        viewUp[1] = 1.0 - upProjection*viewOut[1];
        viewUp[2] =  - upProjection*viewOut[2];

        // Check for validity:
        double upMagnitude = viewUp[0]*viewUp[0] + viewUp[1]*viewUp[1] + viewUp[2]*viewUp[2];

        if (upMagnitude < .0000001) {
            //Second try at making a View Up vector: Use Y axis default  (0,1,0)
            viewUp[0] = -viewOut[1]*viewOut[0];
            viewUp[1] = 1-viewOut[1]*viewOut[1];
            viewUp[2] = -viewOut[1]*viewOut[2];

            // Check for validity:
            upMagnitude = viewUp[0]*viewUp[0] + viewUp[1]*viewUp[1] + viewUp[2]*viewUp[2];

            if (upMagnitude < .0000001) {
                //Final try at making a View Up vector: Use Z axis default  (0,0,1)
                viewUp[0] = -viewOut[2]*viewOut[0];
                viewUp[1] = -viewOut[2]*viewOut[1];
                viewUp[2] = 1-viewOut[2]*viewOut[2];

                // Check for validity:
                upMagnitude = viewUp[0]*viewUp[0] + viewUp[1]*viewUp[1] + viewUp[2]*viewUp[2];

                if (upMagnitude < .0000001) {
                    throw new RuntimeException("Couldn't not find appropriate up vector");
                }
            }
        }

        // normalize the Up Vector
        upMagnitude = Math.sqrt(upMagnitude);
        viewUp[0] = viewUp[0]/upMagnitude;
        viewUp[1] = viewUp[1]/upMagnitude;
        viewUp[2] = viewUp[2]/upMagnitude;

        // Calculate the Right Vector. Use cross product of Out and Up.
        viewRight[0] = -viewOut[1]*viewUp[2] + viewOut[2]*viewUp[1];
        viewRight[1] = -viewOut[2]*viewUp[0] + viewOut[0]*viewUp[2];
        viewRight[2] = -viewOut[0]*viewUp[1] + viewOut[1]*viewUp[0];

        // Plug values into rotation matrix A
        // to create R for getting a space frame vector from body frame
        a11=viewRight[0];
        a21=viewRight[1];
        a31=viewRight[2];

        a12=viewUp[0];
        a22=viewUp[1];
        a32=viewUp[2];

        a13=viewOut[0];
        a23=viewOut[1];
        a33=viewOut[2];
    }

    protected void convertToBodyFrame(IVector v) {
        double tt;
        tt= a11*v.x(0) + a21*v.x(1) + a31*v.x(2);
        tvect.setX(0,tt);
        tt=a12*v.x(0) + a22*v.x(1) + a32*v.x(2);
        tvect.setX(1,tt);
        tt = a13*v.x(0) + a23*v.x(1) + a33*v.x(2);
        tvect.setX(2,tt);
        v.E(tvect);
    }

    //V_space = A*V_body
    protected void convertToSpaceFrame(IVector v) {
        double tt;

        tt = a11*v.x(0) + a12*v.x(1) + a13*v.x(2);
        tvect.setX(0,tt);
        tt= a21*v.x(0) + a22*v.x(1) + a23*v.x(2);
        tvect.setX(1,tt);
        tt = a31*v.x(0) + a32*v.x(1) + a33*v.x(2);
        tvect.setX(2,tt);
        v.E(tvect);
    }

    public double getSigma() {return sigma;}

    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s*s;
    }

    public double getEpsilon() {return epsilon;}

    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4*epsilon;
    }

    public final void setQuadrupolarMomentSquare(double moment){
        Q2=moment;	
    }

    private static final long serialVersionUID = 1L;
    private double sigma , sigma2;
    private double epsilon, epsilon4;
    private double hsdiasq=1.0/Math.sqrt(2);
    private double Q2;
    private NearestImageTransformer nearestImageTransformer;

    double a11, a12, a13, a21,a22,a23,a31,a32,a33;

    private final IVector com1;
    private final IVector com2;
    private final IVector v12;
    private final IVector tvect;
    private final IVector v1;
    private final IVector v2;
    private final IVector az1;
    private final IVector az2;
    private final IVector dr;
    //rotation matrix parameter
    double viewOut[]=new double[3];      // the View or "new Z" vector
    double viewUp[]= new double[3];       // the Up or "new Y" vector
    double viewRight[]=new double[3];    // the Right or "new X" vector

    double worldUp[] = new double[3]; // y axis
}