package etomica.models.oneDHardRods;

/**
 * Class which is used to store which wave vectors we want to compare systems 
 * as harmonic to hard rod potential.
 * 
 * @author cribbin
 *
 */

public class AffectedWaveVectors {
    boolean[] useThisWaveVector;
    int length;
    int numberAffected;
    
    public AffectedWaveVectors(int l){
        length = l;
        useThisWaveVector = new boolean[l];
    }

    /**
     * Argument is a list of the wavevectors that are compared.
     * @param values
     */
    public void setWVs(int[] values){
        numberAffected = 0;
        if(values.length > length){
            throw new IllegalArgumentException("Trying to compare more " +
                    "wavevectors than the system has (in AffectedWaveVectors!");
        }
        //set all modes to used
        for (int i = 0; i < length; i++){
            useThisWaveVector[i] = true;
        }
        //change the modes that are not used to false
        for(int j = 0; j < values.length; j++){
            if(values[j] < 1){
                throw new IllegalArgumentException ("AffectedWaveVectors " +
                        "cannot cope with a negative or zero wavevector");
            }
            useThisWaveVector[values[j]] = false;
            numberAffected++;
        }
    }

    public boolean isUsed(int i){
        return useThisWaveVector[i];
    }

    public boolean[] getWVs(){
        return useThisWaveVector;
    }
    
    public int getNumberOfWVsAffected(){
        return numberAffected;
    }
    
    
//    public static void main(String[] args){
//        int dork = 10;
//        AffectedWaveVectors wvs = new AffectedWaveVectors(dork);
//        
//        System.out.println("1, 6, 5");
//        wvs.setWVs(new int[] {1, 6, 5});
//        for(int j = 0; j < dork; j++){
//            System.out.println("j:  "+j+ "  "+ wvs.getWV(j));
//        }
//        
//        System.out.println("1, 6, 5, 10");
//        wvs.setWVs(new int[] {1, 6, 5, 10});
//        for(int j = 0; j < dork; j++){
//            System.out.println("j:  "+j+ "  "+ wvs.getWV(j));
//        }
//        
//        System.out.println("1, 0");   
//        wvs.setWVs(new int[] {1, 0});
//        for(int j = 0; j < dork; j++){
//            System.out.println("j:  "+j+ "  "+ wvs.getWV(j));
//        }
//        
//        System.out.println("2, -1");
//        wvs.setWVs(new int[] {2, -1});
//        for(int j = 0; j < dork; j++){
//            System.out.println("j:  "+j+ "  "+ wvs.getWV(j));
//        }            
//    }  
}

