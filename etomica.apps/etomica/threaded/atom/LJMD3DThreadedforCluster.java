package etomica.threaded.atom;



import etomica.threaded.IntegratorVelocityVerletThreaded;





public class LJMD3DThreadedforCluster {

    public static void main(String[] args) {
        
        // On-the-fly input - Number of Atoms
        int n;
        n = Integer.parseInt(args[0]);
        
        // On-the-fly input - Number of Threads
        int t;
        t = Integer.parseInt(args[1]);
        
        // On-the-fly input - Number of Timesteps
        int s;
        s = Integer.parseInt(args[2]);
                
        LJMD3DThreaded sim = new LJMD3DThreaded(n, t);
               
       
        sim.activityIntegrate.setMaxSteps(s);
      
        
        // Timer
        double time1 = System.currentTimeMillis();
        
        sim.getController().actionPerformed();
        
        double time2 = System.currentTimeMillis();
        time2 = (time2 - time1)/1000;
        
        
        
        // Output
        System.out.println(t+"_thread(s) - "+n+"_atoms - "+s+"_timesteps");
        

        //System.out.println("force sum time: "+PotentialCalculationForceSumThreaded.forceSumTime/1000+" seconds.");
        
      //  for(int i=0;i<((PotentialMasterListThreaded)sim.getPotentialMaster()).threads.length; i++){
      //  System.out.println("thread"+i+" calculate time: "+((PotentialMasterListThreaded)sim.getPotentialMaster()).threads[i].threadCalculate /1000+" seconds.");
        // }
        
        System.out.println("threads time: "+IntegratorVelocityVerletThreaded.threadtime/1000+" seconds.");
        System.out.println("total time: "+time2+" seconds.");
    }

   
}