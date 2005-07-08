package etomica.utility;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.util.HashMap;
import java.util.LinkedList;

/**
 * Extends the java ObjectInputStream to add fields and methods to enable complete
 * reconstruction of simulation from a serialized form.
 */

/*
 * Created on Jul 6, 2005
 */

public class EtomicaObjectInputStream extends ObjectInputStream {

    public LinkedList atomLists = new LinkedList();
    public HashMap linkerLists = new HashMap();
    
    /**
     * @param out
     * @throws java.io.IOException
     */
    public EtomicaObjectInputStream(InputStream out) throws IOException {
        super(out);
    }

    /**
     * @throws java.io.IOException
     * @throws java.lang.SecurityException
     */
    public EtomicaObjectInputStream() throws IOException, SecurityException {
        super();
    }

}
