package etomica.performance;



import etomica.*;



/**



 * Class for calculating the energy in a procedural way.



 */







public class PotentialBase extends Potential



{

    private double coord[][];

    private double dimensions[];

    private boolean all=false;

    private int index;

    private int N;

  

    /*

    public void PotentialBase(ParentPotential parent, int nAtoms){



        super(parent)

        N=nAtoms;

        coord = new double[N][];

        dimensions = new double[3];



    }

    */

    public  PotentialBase(PotentialGroup parent, int nAtoms){



	    super(parent);

        N=nAtoms;

        coord = new double[N][];

        dimensions = new double[3];



    }

    

    public  Potential set(Atom a){return this;}

    public  Potential set(Atom a1, Atom a2){return this;}

    public  Potential set(SpeciesMaster s){return this;}    

    

    

    public void setArray(Phase phase){

         dimensions[0]=phase.boundary().dimensions().component(0);
         dimensions[1]=phase.boundary().dimensions().component(1);
         dimensions[2]=phase.boundary().dimensions().component(2);

        initialize(phase.makeAtomIterator());   

    }

    

    /**

    * Setting the pointer of coord array to atom.coord.X

    */ 

    

    public void initialize(AtomIterator iterator){

	iterator.reset();

        for(int i=0;i<N;i++){

          coord[i] = iterator.next().coord.position().toArray();

        }

    }



    //another approach, which makes new arrays for coordinates

    public void alternativeInitialize(AtomIterator iterator) {

	iterator.reset();

	for(int i=0; i<N; i++) {

	    Atom atom = iterator.next();

	    double[] X = new double[3];

	    for(int j=0; j<3; j++) {

		X[j] = atom.coord.position().toArray()[j];

	    }

	    coord[i] = X;

	    ((SpaceP.Vector)atom.coord.position()).setArray(X);

	}

    }

    

    public void reinitializeAtomCoord(AtomIterator iterator){

        

       for(int i=0;i<N;i++) iterator.next().coord.position().E(coord[i]);

        

    }

    

    public void reset(IteratorDirective id){

    

        switch(id.atomCount()) {

            case 0: all=true;



                     break;



            case 1:  index=id.atom1().index(); 

                all = false;



                break;



        }



    }



    



    /**



    * Calculate Method for doing energy sum



    */



    public void calculate(IteratorDirective id, PotentialCalculation pc) {

             reset(id);
        double sum=0.0;

        if(all){
            for(int i=0; i<N;i++)
                for(int j=i+1;j<N;j++){
        		    sum+=energy(i,j);



		   // if(sum >= Double.MAX_VALUE) return sum;    
                }
                
           //     System.out.println("sum .. " + sum);
            }



        else {


           
            for(int i=0;i<N;i++){



             if(i!=index) sum+=energy(i,index);



           //  if(sum >= Double.MAX_VALUE) return sum;



            }


            

        }


 
///        ((etomica.PotentialCalculationEnergySumPerformance)pc).setSum(sum);

        



    }//end of calculate



    



    public double energy(int i, int j){



	double[] ci = coord[i];

	double[] cj = coord[j];



	double dx = ci[0] - cj[0];

	double dy = ci[1] - cj[1];

	double dz = ci[2] - cj[2];



	/*

        double dx=coord[i][0]-coord[j][0];  



        double dy=coord[i][1]-coord[j][1];



        double dz=coord[i][2]-coord[j][2];

	*/

        



	//PBC



	//   (could store 0.5*dimensions[i])



	if(dx > 0.5*dimensions[0]) dx -= dimensions[0];



	else if(dx < -0.5*dimensions[0]) dx += dimensions[0];







	if(dy > 0.5*dimensions[1]) dy -= dimensions[1];



	else if(dy < -0.5*dimensions[1]) dy += dimensions[1];







	if(dz > 0.5*dimensions[2]) dz -= dimensions[2];



	else if(dz < -0.5*dimensions[2]) dz += dimensions[2];







        



       // LJ

        return u(dx*dx+dy*dy+dz*dz);

        

    }



    



    //moved this into energy



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

    public double u(double r2) {
          double enrgy; 
	        if(r2 != r2Last) {

            double s2 = sigmaSquared/r2;

            s6 = s2*s2*s2;

            r2Last = r2;

	            }
        enrgy=epsilon4*s6*(s6 - 1.0);
        if(enrgy >Double.MAX_VALUE)
        System.out.println("enrgy " + enrgy);
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



   

    private double sigma, sigmaSquared;

    private double epsilon;

    private double epsilon4, epsilon48, epsilon624;

    private static final double _168div624 = 168./624.;

    private double uInt, rCLast;  

    private double r2Last = -1.0;

    private double s6;



}



