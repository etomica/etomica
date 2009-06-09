package etomica.virial.cluster2.nauty;

import java.io.IOException;
import java.io.Reader;
import java.util.Map;

public interface ProcessWrapper {

  public Map<String, String> getEnvironment();

  public ProcessInfo getProcessInfo();

  public Reader getProcessOutput();

  public void run() throws IOException;
}