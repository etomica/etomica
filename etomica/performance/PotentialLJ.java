package etomica.performance;
import etomica.*;

public class PotentialLJ extends PotentialBase{
    private double sigma, sigmaSquared;
    private double epsilon;
    private double epsilon4, epsilon48, epsilon624;
    private static final double _168div624 = 168./624.;
    private double uInt, rCLast;  
    private double r2Last = -1.0;
    private double s6;
    private Space3D.Vector v = new Space3D.Vector();
    double a[]=new double[3];
    SpaceP.BoundaryPeriodicSquare boundary;
    
    public PotentialLJ(PotentialGroup parent, int nAtoms){
        super(parent,nAtoms);
       // v=(SpaceP.Vector)parentSimulation().space().makeVector();
    }
    
    public void setArray(Phase phase1) {
        super.setArray(phase1);
        boundary = (SpaceP.BoundaryPeriodicSquare)phase1.boundary();
    }
    
    public final double energy(int l, int m){
        double[] cl = coord[l];
        double[] cm = coord[m];
        double dx = cl[0] - cm[0];
        double dy = cl[1] - cm[1];
        double dz = cl[2] - cm[2];
       
     //   double dx = coord[l][0] - coord[m][0];
     //   double dy = coord[l][1] - coord[m][1];
     //   double dz = coord[l][2] - coord[m][2];
        
   /*     while(dx > d02x)dx-= dmx;
        while(dx <_d02x)dx+= dmx;
        while(dy > d02y)dy-=dmy;
        while(dy <_d02y)dy+=dmy;
        while(dz > d02z)dz-=dmz;
        while(dz <_d02z)dz+=dmz;
       
        return u(dx*dx+dy*dy+dz*dz);
     */    
        
        //double a[]=new double[3];
        
     //   a[0]=dx;a[1]=dy;a[2]=dz;
    //    v.E(a);
       // Space.Vector v = parentSimulation().space().makeVector();
       //v.setComponent(0,dx);
       
       
        v.E(dx,dy,dz);
        
        //v.E(a);
        //return u(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
        
   ///       boundary.nearestImage(a);
         boundary.nearestImage(v);
      //  return u(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
        
     //  return u(dx*dx+dy*dy+dz*dz);
       return u(v.squared());

    }
    
    public double u(double r2) {
          //double enrgy; 
//	       if(r2 != r2Last) {
    
            double s2 = sigmaSquared/r2;

            s6 = s2*s2*s2;

   //         r2Last = r2;

	//            }
       //enrgy=epsilon4*s6*(s6 - 1.0);
      //  if((all)&&(s6 > 1.0)) System.out.println("panga in u" + enrgy );
        return epsilon4*s6*(s6 - 1.0);

    }
    /**
     * Mutator method for Lennard-Jones size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */

    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
        rCLast = 0.0;
    }

   public final void setEpsilon(double eps) {

        uInt *= eps/epsilon;
        epsilon = eps;
        epsilon4 = eps*4.0;
        epsilon48 = eps*48.0;
        epsilon624 = eps*624.0;
    }

    
    
    
}
