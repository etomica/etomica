package etomica.virial.cluster.graphics;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import etomica.virial.cluster.*;


/**
 * @author skkwak
 * 
 *  simple method to show cluster
 * 	it takes number of points and bond lists
 */

public class ShowCluster extends JFrame {
	private int nPoints;
    private JPanel compositePanel;
    private JPanel controlPanel;
    private JButton forwardButton;
    private JButton backwardButton;
    private JButton resetButton;
    private JButton toggleButton;
    private ClusterPanel clusterPanel;
    private ControlListener cl;
	private double[] x;
	private double[] y;
    private int[][] bonds;
	private java.util.LinkedList testList = new java.util.LinkedList();
	
	public ShowCluster(int nPoints, java.util.LinkedList testList){
		super( " Cluster : Total Number of Clusters = "+Integer.toString(testList.size()));
		setSize( 450, 500);
		Container contentPane = getContentPane();
		this.nPoints = nPoints;
		this.testList = testList;
		cl = new ControlListener();
		cl.setTestList(testList);
		
        contentPane.add(getCompositePanel());
		show();
	}
	
	public JPanel getCompositePanel(){
		compositePanel = new JPanel();
		controlPanel = new JPanel();
		controlPanel.add(resetButton = new JButton("reset"));
			resetButton.addActionListener(cl);
		controlPanel.add(backwardButton = new JButton("Backward"));
			backwardButton.addActionListener(cl);
		controlPanel.add(forwardButton = new JButton("Forward"));
			forwardButton.addActionListener(cl);
		controlPanel.add(toggleButton = new JButton("Toggle"));
			toggleButton.addActionListener(cl);
		clusterPanel = new ClusterPanel();
		compositePanel.add(clusterPanel);		
		compositePanel.add(controlPanel);
		return compositePanel;		
	}
		
	public class ClusterPanel extends Plot {
    	private double max=10.0;
    	private double m = 0;
    	private double k = 0;
    	private boolean needToClear = true;
		private boolean clearPoints = false;

		public ClusterPanel(){
	     setSize(400,400);
	     setBackground("white");
		}
		public void paint(){
	      
	      if(needToClear){
	      clear();
	      needToClear = false;
    	  double size=0.5;
    	  double r = 5.0;
    	  double r1 =7.0;
    	  double rotateAngle = 2.0*Math.PI/nPoints;
    	  double angle = -Math.PI/2.0+(rotateAngle/2.0);
    	  x = new double[nPoints];
    	  y = new double[nPoints];
	      setColor("white");
    	  drawAxes(-max,max,-max,max);
    	  setColor("black");
	    	  for(int i=0; i<nPoints;i++){
	    	  		x[i]=-r*Math.cos(angle+(i*rotateAngle));
	    	  		y[i]=r*Math.sin(angle+(i*rotateAngle));
	    	  	if(i==0 || i==(nPoints-1)){ //setColor("red");	    	  		 
		    	  if(!clearPoints){	circle(x[i]-size,x[i]+size,y[i]-size,y[i]+size);
		    	  					plotStringCenter(Integer.toString(i),-r1*Math.cos(angle+(i*rotateAngle)),r1*Math.sin(angle+(i*rotateAngle)));
		    	  }
	    	  	} else { //setColor("blue");
		    	  if(!clearPoints){	floodCircle(x[i]-size,x[i]+size,y[i]-size,y[i]+size);
		    	  	plotStringCenter(Integer.toString(i),-r1*Math.cos(angle+(i*rotateAngle)),r1*Math.sin(angle+(i*rotateAngle)));
		    	  }
	    	  	}//end of if
	    	  }//end of for loop
	      }//end of outer if need to clear
	     
	     setColor("black");	
	     plotOverlay();
	     bonds = (int[][])cl.getBondGroup();
		 plotStringCenter("GN = "+cl.testList.indexOf(cl.getBondGroup())+
		                  "   NoB = "+bonds.length,0,max);	     
	     for(int i=0;i<bonds.length;i++){	 
	 	     plotLine(x[bonds[i][0]],y[bonds[i][0]] , x[bonds[i][1]],y[bonds[i][1]]);
             plotStringCenter("["+bonds[i][0]+" "+bonds[i][1]+"]",-(max+1),((max-1)-i));
	     }
		}//end of paint()	
		
		public void setNeedToClearPoints(boolean clearPoints){
			needToClear = true;
			this.clearPoints = clearPoints;
		}
		public boolean getNeedToClearPoints(){
			return clearPoints;
		}
	}// end of ClusterPanel
		
	public class ControlListener implements ActionListener {
		public java.util.LinkedList testList = new java.util.LinkedList();
        private Object bondGroup;
		private int index = 0;

        public void setTestList(java.util.LinkedList testList){
			this.testList = testList;        	
			setBondGroup(testList.getFirst());
            setIndex(testList.indexOf(testList.getFirst()));
        }
        public void setBondGroup(Object bonds){
        	bondGroup = bonds;
        }
        public Object getBondGroup(){
        	return bondGroup;
        }
        public void setIndex(int index){
        	this.index = index;
        }
        public int getIndex(){
        	return index;
        }
        
		public void actionPerformed(ActionEvent event){
			if(event.getSource()==resetButton){
				
				setBondGroup(testList.getFirst());
                setIndex(testList.indexOf(testList.getFirst()));
		        clusterPanel.repaint();

			} else if(event.getSource()==backwardButton){

		        int n = getIndex()-1;
		        if(n>=0){ setBondGroup(testList.get(n));
		        	      setIndex(n);
		        	      clusterPanel.repaint();
		        } else {System.out.println(" First Cluster !!!! ");}			

			} else if(event.getSource()==forwardButton){

		        int n = getIndex()+1;
		        if(n<=testList.indexOf(testList.getLast())){ setBondGroup(testList.get(n));
		        	      setIndex(n);
		        	      clusterPanel.repaint();
		        } else {System.out.println(" Last Cluster !!!! ");}			

			} else if(event.getSource()==toggleButton){
				 if(clusterPanel.getNeedToClearPoints()){
				  clusterPanel.setNeedToClearPoints(false);
				 } else {clusterPanel.setNeedToClearPoints(true);}
				 clusterPanel.repaint();
			}
		}
	}

/*
	public class TotalClusterPanel extends Plot{
		private java.util.LinkedList testList = new java.util.LinkedList();
		private java.util.Iterator iterator;
    	private double max=10.0;
    	private double m = 0;
    	private double k = 0;
		public TotalClusterPanel(){
	     setSize(400,400);
	     setBackground("white");
		}
		public void paint(){
	      if(needToClear){
	      clear();
	      needToClear = false;
	      
    	  double size=0.5;
    	  double r = 2.0;
    	  double r1 =3.0;
    	  double rotateAngle = 2.0*Math.PI/nPoints;
    	  double angle = -Math.PI/2.0+(rotateAngle/2.0);
    	  x = new double[nPoints];
    	  y = new double[nPoints];
	      setColor("white");
    	  drawAxes(-max,max,-max,max);
    	  setColor("black");
    	  while(iterator.hasNext()){
    	  	  bonds = (int[][])iterator.next();
	    	  for(int i=0; i<nPoints;i++){
	    	  		x[i]=-r*Math.cos(angle+(i*rotateAngle));
	    	  		y[i]=r*Math.sin(angle+(i*rotateAngle));
	    	  	if(i==0 || i==(nPoints-1)){ //setColor("red");	    	  		 
		    	  	circle(x[i]-size,x[i]+size,y[i]-size,y[i]+size);
//		    	  	plotStringCenter(Integer.toString(i),-r1*Math.cos(angle+(i*rotateAngle)),r1*Math.sin(angle+(i*rotateAngle)));
	    	  	} else { //setColor("blue");
		    	  	floodCircle(x[i]-size,x[i]+size,y[i]-size,y[i]+size);
//		    	  	plotStringCenter(Integer.toString(i),-r1*Math.cos(angle+(i*rotateAngle)),r1*Math.sin(angle+(i*rotateAngle)));
	    	  	}
			     for(int j=0;i<bonds.length;j++){	 
			 	     plotLine(x[bonds[j][0]],y[bonds[i][0]] , x[bonds[j][1]],y[bonds[j][1]]);
//		             plotStringCenter("["+bonds[j][0]+" "+bonds[j][1]+"]",-(max+1),((max-1)-j));
			     }//end of inner For loop	    	  	
	    	  }//end of outer For loop
    	  }// end of while

	      }//end of outer if
	     
		}//end of paint()	

        public void setTestList(java.util.LinkedList testList){
			this.testList = testList;  
			iterator = testList.iterator();      	
        }
	}// end of TotalClusterPanel
*/
	
	public static void main(String[] args) {
//		ShowCluster sc = new ShowCluster(4);
//		sc.addWindowListener(new WindowAdapter(){
//			public void windowClosing( WindowEvent e){
//				System.exit(0);
//			}});
	}
}
