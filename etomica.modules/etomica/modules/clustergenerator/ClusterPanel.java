/*
 * Created on Dec 9, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.clustergenerator;

import etomica.virial.cluster.ClusterDiagram;

//import etomica.virial.simulations.GenCluster;


public class ClusterPanel extends Plot {
    private double[] x;
    private double[] y;
    private boolean drawPoints;
    private boolean drawNumbersWhite, drawNumbersBlack;
    private final boolean fBonds, eBonds;
    private ClusterDiagram cluster;

    public ClusterPanel(ClusterDiagram aCluster,boolean drawPoints, boolean drawNumbersWhite, boolean drawNumbersWBlack, 
            boolean fDrawBonds, boolean eDrawBonds) {
        super(null,400,400);
        cluster = aCluster;
        x = new double[cluster.mNumBody];
        y = new double[cluster.mNumBody];
        this.drawPoints = drawPoints;
        this.drawNumbersWhite = drawNumbersWhite;
        this.drawNumbersBlack = drawNumbersWBlack;
        fBonds = fDrawBonds;
        eBonds = eDrawBonds;
    }
    
    public ClusterDiagram getCluster(){
    	return cluster;
    }

    public void paint(){
        int nPoints = cluster.mNumBody;
        clear();
        double size=1.3;
        double max=13.0;
        double r = 8.0;
        double rotateAngle = 2.0*Math.PI/nPoints;
        double angle = -Math.PI/2.0 - rotateAngle*((cluster.getNumRootPoints()-1)*0.5);
        setWindow(-max,max,-max,max);
        for(int i=0; i<nPoints;i++){
            x[i]=-r*Math.cos(angle+(i*rotateAngle));y[i]=r*Math.sin(angle+(i*rotateAngle));
        }

        for (int i=0; i<nPoints; i++) {
            int lastBond = i;
            int[] iConnections = cluster.mConnections[i];
            for (int j=0; j<nPoints-1; j++) {
                if (iConnections[j] > i) {
                    if (eBonds) {
                        setColor("red");
                        for (int k=lastBond+1; k<iConnections[j]; k++) {
                            if (!cluster.isRootPoint(i) || !cluster.isRootPoint(k)) {
                                plotLine(x[i],y[i], x[k],y[k]);
                            }
                        }
                        setColor("black");
                    }
                    if (fBonds) {
                        plotLine(x[i],y[i] , x[iConnections[j]],y[iConnections[j]]);
                    }
                    lastBond = iConnections[j];
                }
                else if ((lastBond>i || iConnections[j] == -1) && eBonds) {
                    setColor("red");
                    for (int k=lastBond+1; k<nPoints; k++) {
                        if (!cluster.isRootPoint(i) || !cluster.isRootPoint(k)) {
                            plotLine(x[i],y[i], x[k],y[k]);
                        }
                    }
                    setColor("black");
                }
                if (iConnections[j] == -1) break;
            }
        }

        setColor("black");
        for(int i=0; i<nPoints;i++){  
            if (drawPoints) {
                if(i < cluster.getNumRootPoints()) {
                    setColor("white");
                    floodInside(x[i]-size,x[i]+size,y[i]-size,y[i]+size);
                    setColor("black");
               	    circle(x[i]-size,x[i]+size,y[i]-size,y[i]+size);
                } else {
                    floodCircle(x[i]-size,x[i]+size,y[i]-size,y[i]+size);
                }
            }
            if(i < cluster.getNumRootPoints() && drawNumbersWhite) {
                plotStringCenter(Integer.toString(i),(-r*1.4)*Math.cos(angle+(i*rotateAngle))-0.2,(r*1.4)*Math.sin(angle+(i*rotateAngle))-1.0);               
            } else if(i >= cluster.getNumRootPoints() && drawNumbersBlack){
                plotStringCenter(Integer.toString(i),(-r*1.4)*Math.cos(angle+(i*rotateAngle))-0.2,(r*1.4)*Math.sin(angle+(i*rotateAngle))-1.0);                          
            }
        }
	}	
}