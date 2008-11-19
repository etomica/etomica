package etomica.models.oneDHardRods;

/**
 * Class which is used to store which wave vectors we want to compare systems 
 * as harmonic to hard rod potential.
 * 
 * @author cribbin
 *
 */

public class AffectedWaveVectors {
    boolean[] useme;
    int length;
    
    public AffectedWaveVectors(int l){
        length = l;
        useme = new boolean[l];
    }

    /**
     * Argument is a list of the wavevectors that are affected.
     * @param values
     */
    public void setWVs(int[] values){
        if(values.length < length){
            throw new IllegalArgumentException("Trying to compare more " +
                    "wavevectors than the system has (in AffectedWaveVectors!");
        }
        for (int i = 0; i < length; i++){
            useme[i] = false;
        }
        for(int j = 0; j < values.length; j++){
        }
    }

    public boolean getWV(int i){
        return useme[i];
    }

    public boolean[] getWVs(){
        return useme;
    }
}
