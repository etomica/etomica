package etomica.performance;
import etomica.*;
import etomica.units.Dimension;
public class PotentialHardAssociationCone extends PotentialBase{
    private double wellCutoff, wellCutoffSquared;
    private double sigma, sigmaSquared;
    private double epsilon, epsilon4, wellEpsilon;
    private double cutoffLJSquared, cutoffFactor;
    private Space.Vector e1;
    private Space.Vector e2;
    private double theta, ec2,s2,s6;
    private double v[];
    int i,j;
   
    public PotentialHardAssociationCone(PotentialGroup parent, int nAtoms){
        super(parent,nAtoms);
        v=new double[3];
        e1 = parentSimulation().space().makeVector();
        e2 = parentSimulation().space().makeVector();

        setSigma(Default.ATOM_SIZE);
        setEpsilon(Default.POTENTIAL_WELL);
        
        setCutoffFactorLJ(Default.POTENTIAL_CUTOFF_FACTOR);
        setWellCutoff(getSigma());
        setWellEpsilon(8.0*getEpsilon());
        setTheta(etomica.units.Degree.UNIT.toSim(27.0));
        
    }  
        
      
    public double energy(int l, int m){
        i=l;j=m;
        double dx = coord[l][0] - coord[m][0];
        double dy = coord[l][1] - coord[m][1];
    	double dz = coord[l][2] - coord[m][2];
    	
        v[0]=dx;v[1]=dy;v[2]=dz;
    	
        while(dx > d02x)dx-= dmx;
        while(dx <_d02x)dx+= dmx;
        while(dy > d02y)dy-=dmy;
        while(dy <_d02y)dy+=dmy;
        while(dz > d02z)dz-=dmz;
        while(dz <_d02z)dz+=dmz;
        return u(dx*dx+dy*dy+dz*dz);
        //return u(l,m);
    }
    
    //public double u(double dx, double dy, double dz,int i, int j) {
    //public double u(int i, int j) {
    public double u(double r2) {
       // r2=dx*dx+dy*dy+dz*dz;
      double eTot = 0.0;
        
       
        
        
        if(r2 > cutoffLJSquared) {
            eTot = 0.0;
        }
        else {
            
            s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            eTot = epsilon4*s6*(s6 - 1.0);
        }
        
        
                
        if (r2 < wellCutoffSquared) {
            e1.E(0.);
            
            e1.setComponent(0,1.0);
            
            e1.setComponent(0,e1.component(0)*rotmatrix[i][0][0]);
            e1.setComponent(1,e1.component(1)*rotmatrix[i][1][0]);
            e1.setComponent(1,e1.component(2)*rotmatrix[i][2][0]);
            
           
           double er1 = e1.dot(v);
                       
            if ( er1 > 0.0 && er1*er1 > ec2*r2) {
                e2.E(0.);
                e2.setComponent(0,1.0);
                
                e2.setComponent(0,e2.component(0)*rotmatrix[j][0][0]);
                e2.setComponent(1,e2.component(1)*rotmatrix[j][1][0]);
                e2.setComponent(1,e2.component(2)*rotmatrix[j][2][0]);
                
                double er2 = e2.dot(v);
                if(er2 < 0.0 && er2*er2 > ec2*r2) eTot -= wellEpsilon;
            }
        }
        
        
        return eTot;  
         
         
    }
    
      /**
     * Accessor method for Lennard-Jones size parameter
     */
    public double getSigma() {return sigma;}
    /**
     * Accessor method for Lennard-Jones size parameter
     */
    public void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
        setCutoffFactorLJ(cutoffFactor);
    }
    public static final Dimension getSigmaDimension() {return Dimension.LENGTH;}

    /**
    * Accessor method for Lennard-Jones cutoff distance; divided by sigma
    * @return cutoff distance, divided by size parameter (sigma)
    */
    public double getCutoffFactorLJ() {return cutoffFactor;}
    /**
     * Accessor method for Lennard-Jones cutoff distance; divided by sigma
     * @param rc cutoff distance, divided by size parameter (sigma)
     */
    public void setCutoffFactorLJ(double rc) {  
        cutoffFactor = rc;
        double cutoffLJ = sigma*cutoffFactor;
        cutoffLJSquared = cutoffLJ*cutoffLJ;
       // calculateLRC();
    }
    public static final Dimension getCutoffFactorLJDimension() {return Dimension.NULL;}
   
    /**
    * Accessor method for Lennard-Jones energy parameter
    */ 
    public double getEpsilon() {return epsilon;}
    /**
    * Accessor method for depth of well
    */
    public void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4.0 * eps;
    }
    public static final Dimension getEpsilonDimension() {return Dimension.ENERGY;}
    
    /**
    * Accessor method for attractive-well diameter.
    */
    public double getWellCutoff() {return wellCutoff;}
    /**
    * Accessor method for attractive-well diameter.
    */
    public void setWellCutoff(double wcut) {
        wellCutoff = wcut;
        wellCutoffSquared = wcut*wcut;
    }
          
    public static final Dimension getWellCutoffDimension() {return Dimension.LENGTH;}
    
    /**
    * Accessor method for attractive-well depth parameter.
    */
    public double getWellEpsilon() {return wellEpsilon;}
    /**
    * Accessor method for attractive-well depth parameter.
    */
    public void setWellEpsilon(double weps) {wellEpsilon = weps;}
          
    public static final Dimension getWellEpsilonDimension() {return Dimension.ENERGY;}
    
    /**
     * Accessor method for angle describing width of cone.
     */
    public double getTheta() {return theta;}
    
    /**
     * Accessor method for angle (in radians) describing width of cone.
     */
    public void setTheta(double t) {
        theta = t;
        ec2    = Math.cos(theta);
        ec2   = ec2*ec2;
    }
    public Dimension getThetaDimension() {return Dimension.ANGLE;}
    
    
}