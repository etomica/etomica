package etomica.models.oneDHardRods;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.IAction;
import etomica.util.Histogram;

public class WriteHistograms implements IAction{
    
    Histogram histogram; 
//    FileWriter writer;
    String filename;
    
    public WriteHistograms(){
        
    }
    
    
    public WriteHistograms(String fn){
        filename = fn;
    }
    
    
    public void actionPerformed(){
        
        try{
            FileWriter writer = new FileWriter(filename);
            
            int nBins = histogram.getNBins();
            
            for(int i = 0; i < nBins; i++){
                writer.write(Double.toString(histogram.xValues()[i]));
            }
            
            writer.write("\n");
            
            for(int i = 0; i < nBins; i++){
                writer.write(Double.toString(histogram.getHistogram()[i]));
            }
            
        } catch (IOException e) {
            throw new RuntimeException("Oops, failed to write data " + e);
        }
    }

    public void setFilename(String fn){
        filename = fn;

    }
    
    public void setHistogram(Histogram h){
        this.histogram = h;
    }
}
