/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.IAction;
import etomica.data.histogram.Histogram;

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
                writer.write(Double.toString(histogram.xValues()[i]) + " " + Double.toString(histogram.getHistogram()[i]) +"\n");
            }
            
            writer.close();
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
