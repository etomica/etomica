package etomica.gui;

import java.awt.*;
import java.awt.image.*;
import java.net.URL;
import javax.swing.*;
import java.io.Serializable;
import java.io.InputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;
import java.io.IOException;

/** A simple implementation of the <code>Icon</code> interface.
Similar to Sun's <a href="com.sun.java.swing.ImageIcon"><code>ImageIcon</code></a>.
*/
public class ImageIcon
	implements Icon, java.io.Serializable
{	
	/**
	 * Creates an uninitialized image icon.
	 */
	public ImageIcon()
	{
	}
	
	/**
	 * Creates an image icon from the specified URL. The image will
	 * be preloaded by using MediaTracker to monitor the loaded state
	 * of the image.
	 */
	public ImageIcon (URL location)
	{
		imageLocation = location;
		image = Toolkit.getDefaultToolkit().getImage(location);
		loadImage(image);
	}
	
	/**
	 * Creates an image icon from the specified Image. The image will
	 * be preloaded by using MediaTracker to monitor the loaded state
	 * of the image.
	 */
	public ImageIcon (Image i)
	{
	    image = i;
	    loadImage(image);
	}
	
	/**
	 * Creates an image icon from the specified InputStream. The image will
	 * be preloaded by using MediaTracker to monitor the loaded state
	 * of the image.  The InputStream is closed when the read is done.
	 */
	public ImageIcon (InputStream inStream)
	    throws IOException
	{
	    int                   numRead;
	    byte[]                bytes;
	    ByteArrayOutputStream outStream;
	    
	    bytes     = new byte[1024];
	    outStream = new ByteArrayOutputStream();
	    
	    while((numRead = inStream.read(bytes)) > 0)
	    {
	        outStream.write(bytes, 0, numRead);
	    }
	    
	    inStream.close();
	    outStream.flush();
	    
		image = Toolkit.getDefaultToolkit().createImage(outStream.toByteArray());
	    loadImage(image);
	    outStream.close();
	}
	
	//
	// Properties
	//
	
	/**
	 * Returns the location of the Image displayed by this icon.
	 */
	public URL getImageLocation ()
	{
		return imageLocation;
	}
	
	/**
	 * Set the location to displayed by this icon. The image will
	 * be preloaded by using MediaTracker to monitor the loading state
	 * of the image.
	 */
	public void setImageLocation (URL location)
	{
		imageLocation = location;
		image = Toolkit.getDefaultToolkit().getImage(location);
		loadImage(image);
	}
	
	/** 
	 * Set the image observer for the image.  Set this
	 * property if the ImageIcon contains an animated GIF.
	 * For example:
	 * <pre>
	 *     icon = new ImageIcon(...)
	 *     button.setImageLocation(URL);
	 *     icon.setImageObserver(button);
	 * </pre>
	 */
	public void setImageObserver(ImageObserver observer)
	{
		imageObserver = observer;
	}
	
	/**
	 *  Return the image observer for the image 
	 */
	public ImageObserver getImageObserver()
	{
		return imageObserver;
	}
	
	//
	// Icon interface implementation
	//
	
	/**
	 * Paints the Icon
	 */
	public synchronized void paintIcon(Component c, Graphics g, int x, int y)
	{
		if (image != null)
			if (imageObserver == null)
				g.drawImage(image, x, y, c);
			else
				g.drawImage(image, x, y, imageObserver);
	}
	
	/**
	 * Get the width of the Icon
	 */
	public int getIconWidth()
	{
		return width;
	}
	
	/**
	 * Get the height of the Icon
	 */
	public int getIconHeight()
	{
		return height;
	}
	
	//
	// Implementation
	//
	
	/**
	 * Wait for the image to load
	 */
	protected void loadImage(Image image)
	{
		synchronized(tracker)
		{
			tracker.addImage(image, 0);
			
			try
			{
				tracker.waitForID(0, 5000);
			}
			catch (InterruptedException e)
			{
				System.out.println("INTERRUPTED while loading Image");
			}
			
			loadStatus = tracker.statusID(0, false);
			
			tracker.removeImage(image, 0);
			
			width = image.getWidth(imageObserver);
			height = image.getHeight(imageObserver);
		}
	}
	
	/**
	 * Returns the status of the image loading operation.
	 * @return the loading status as defined by java.awt.MediaTracker.
	 * @see java.awt.MediaTracker#ABORTED
	 * @see java.awt.MediaTracker#ERRORED
	 * @see java.awt.MediaTracker#COMPLETE
	 */
	public int getImageLoadStatus()
	{
	    return loadStatus;
	}
	
	/** Gets the <code>Image</code> of this <code>Icon</code>.
	@return the <code>Image</code>
	*/
	public Image getImage()
	{
	    return image;
	}
	
	/** The icon's image.
	@see #getImage
	*/
	protected transient Image image;
	/** The status of image loading.
	@see #getImageLoadStatus
	*/
	protected transient int loadStatus = 0;
	/** The observer of the image.
	@see #getImageObserver
	*/
	protected ImageObserver imageObserver;
	/** The location of the image data.
	@see #getImageLocation
	*/
	protected URL imageLocation = null;
	
	/** Internal use.
	*/
	protected final static Component component = new Component() {};
	/** The tracks image loading.
	*/
	protected final static MediaTracker tracker = new MediaTracker(component);
	
	/** The icon width.
	@see #getIconWidth
	*/
	protected int width = -1;
	/** The icon height.
	@see #getIconHeight
	*/
	protected int height = -1;
}

