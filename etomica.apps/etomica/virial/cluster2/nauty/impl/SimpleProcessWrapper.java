package etomica.virial.cluster2.nauty.impl;

import java.io.*;
import java.util.*;

import etomica.virial.cluster2.nauty.ProcessInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;

public class SimpleProcessWrapper implements ProcessWrapper {

  private Reader output;
  private ProcessBuilder builder;
  private Map<String, String> environment;
  private ProcessInfo processInfo;

  public SimpleProcessWrapper(ProcessInfo info) {

    processInfo = info;
    builder = new ProcessBuilder(info.getCommand());
    environment = builder.environment();
  }

  public Map<String, String> getEnvironment() {

    return environment;
  }

  public Reader getProcessOutput() {

    return output;
  }

// @Override
// public List<String> getStdout() {
//
// return stdout;
// }
  public ProcessInfo getProcessInfo() {

    return processInfo;
  }

  /**
   * Blocking execution of the command.
   * 
   * @throws IOException
   */
  public void run() throws IOException {

    builder.redirectErrorStream(processInfo.redirectErrorStream());
    builder.directory(new File(processInfo.getPath()));
    Process process = builder.start();
    output = new InputStreamReader(process.getInputStream());
  }
}