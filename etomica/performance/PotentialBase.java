package etomica.performance;
import etomica.*;
/**
 * Class for calculating the energy in a procedural way.
 */
public abstract class PotentialBase extends Potential
{

    public double coord[][];
    public double rotmatrix[][][];
    //private double dimensions[];
    public boolean all=false;
    public int index;
    public int N;
    public double invdimx,invdimy,invdimz;
    public double d02x,_d02x, d02y,_d02y,d02z,_d02z;
    public double dmx,dmy,dmz;
    public Phase phase;

    public  PotentialBase(PotentialGroup parent, int nAtoms){
	    super(parent);
        N=nAtoms;
        coord = new double[N][3];
       // dimensions = new double[3];
    }
    
    public  Potential set(SpeciesMaster s){return this;} 
    public  Potential set(Atom[] atoms) {return this;}
    public  void setSpecies(Species[] s) {}
    public  Species[] getSpecies() {return null;}
    public abstract double energy(int l, int m);
    public void setArray(Phase phase1){
        phase=phase1;
        dmx=phase.boundary().dimensions().x(0);
        dmy=phase.boundary().dimensions().x(1);
        dmz=phase.boundary().dimensions().x(2);
        invdimx=1.0/dmx;
        invdimy=1.0/dmy;
        invdimz=1.0/dmz;
        d02x=0.5*dmx;_d02x=-d02x;
        d02y=0.5*dmy;_d02y=-d02y;
        d02z=0.5*dmz;_d02z=-d02z;
        
        initialize(phase.makeAtomIterator());   
        //alternativeInitialize(phase.makeAtomIterator());
    }  
    public void setRotationMatrix(Phase phase){
        AtomIterator iterator = phase.makeAtomIterator();
        iterator.reset();
        rotmatrix=new double[N][3][3];
        
        for(int i=0;i<N;i++){
            Atom a = iterator.next();
            //System.out.println((Space.Coordinate.AngularP)a.coord);
         // System.out.println(" checking ..."+((Space.Coordinate.Angular)a.coord).orientation().getRotMatrix());
            
//           rotmatrix[i]=((Space.Coordinate.AngularP)a.coord).orientation().getRotMatrix();
        }
    }
    /**
    * Setting the pointer of coord array to atom.coord.X
    */ 
    public void initialize(AtomIterator iterator){
        iterator.reset();
        for(int i=0;i<N;i++){
           //System.out.println("array ...." +iterator.next().coord.position().toArray())
            //Atom atom = iterator.next();
            
           coord[i] = iterator.next().coord.position().toArray();
           //coord[i] =atom.coord.position().toArray();
             //  System.out.println(" index " +atom.index() + " i " + i);
          // System.out.println(" checking array initialization..." + coord[i]);
        }
    }
    
    /**
    * Method to put the arrays next to each other. Not sure.
    */
    public void alternativeInitialize(AtomIterator iterator) {
        iterator.reset();
        for(int i=0; i<N; i++) {
            Atom atom = iterator.next();
            double[] X=atom.coord.position().toArray();
            
            for(int j=0; j<3; j++) {
        	    coord[i][j] = X[j];
            }
            //System.out.println(" coord [i] " +coord[i]);
            ((SpaceP.Vector)atom.coord.position()).setArray(coord[i]);
            

	    }

    }

    public void reinitializeAtomCoord(AtomIterator iterator){

       for(int i=0;i<N;i++) iterator.next().coord.position().E(coord[i]);

    }

    public void reset(IteratorDirective id){
        switch(id.atomCount()) {
            case 0: all=true;
                     break;
            case 1:  index=id.atom1().node.index(); 
                all = false;
                break;
        }
    }

    /**
    * Method for calculating energy sum
    */
    public void calculate(IteratorDirective id, PotentialCalculation pc0) {
        reset(id);
        double sum=0.0;
        
        PotentialCalculationEnergySumPerformance pc = (PotentialCalculationEnergySumPerformance)pc0;
        PotentialCalculationEnergySumPerformance.MyAtomPairAction action =
            (PotentialCalculationEnergySumPerformance.MyAtomPairAction)pc.getAtomPairCalculation(this);
        if(all){
            
            for(int i=0; i<N;i++){
                for(int j=i+1;j<N;j++){
        		    sum+=energy(i,j);
                   // action.action(i,j);
		   // if(sum >= Double.MAX_VALUE) return sum;    
                }
            }
               //if(sum >0)
               //System.out.println("sum .. " + sum);
            }

        else {
           // for(int i=0;i<N;i++){
            for(int i=N-1;i>index;i--){
                //action.action(i,index);
                sum+=energy(i,index);
           //  if(sum >= Double.MAX_VALUE) return sum;
            }
            for(int i=index-1;i>=0;i--){
               // action.action(i,index);
               sum+=energy(i,index);
           //  if(sum >= Double.MAX_VALUE) return sum;
            }
        }
        ((etomica.PotentialCalculationEnergySumPerformance)pc).setSum(sum);
    }//end of calculate

 
    
    
    
    public double floor(double val){
     
        if(val >0) return (double)((int)val);
        else return (double)((int)(val-1.0));
        
    }
    
    public double ceil(double val){
            if(val < 0) return (double)(int)val;
            else 
            if((val!=(double)((int)val)))
            return (double)((int)(val+1.0));
            else return  val;
            
       
    }
 /*
    // NOT the best way of doing it.
    public void nearestImage(double[] dr){
        dr[0]=dr[0]/dimensions[0];
        dr[1]=dr[1]/dimensions[1];
        dr[2]=dr[2]/dimensions[2];
        if(dr[0]> 0.5)dr[0] -=1.0; 
        else if(dr[0] <-0.5)dr[0] +=1.0;
        if(dr[1]>1.0) dr[1] -=1.0; 
        else if (dr[1] <-0.5)dr[1] +=1.0;
        if(dr[2]>0.5)dr[2]-=1.0;
        else if (dr[2]<=0.5)dr[2]+=1.0;

        dr[0]*=dimensions[0];
        dr[1]*=dimensions[1];
        dr[3]*=dimensions[2];
     }

*/
   public void setCoordinate(double r[][]){
        coord=r;   

   }
   public double[][] getCoordinate(){

        return coord;
    }
    /**
     * The energy u.  No truncation is applied here; 

     * instead it is applied in the energy(AtomPair) method of Potential2SoftSpherical.

     */




}



