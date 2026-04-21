package etomica.GasMOP;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

class FileData {
    private String fileName;
    private String value1;
    private double value2;
    private double value3;
    private double value5;
    private String value4;

    public FileData(String fileName, String value1, double value2, double value3,  String value4) {
        this.fileName = fileName;
        this.value1 = value1;
        this.value2 = value2;
        this.value3 = value3;
        this.value4 = value4;
    }

    public String getFileName() {
        return fileName;
    }

    public String getValue1() {return value1;}

    public double getValue2() {
        return value2;
    }
    public double getValue3() {
        return value3;
    }

    @Override
    public String toString() {
       // return String.valueOf("value 1" +  value1+"value 2 " + value2 + " value 3 " + value3 ) ;
        return String.valueOf(value2) ;
    }
    public static void main(String[] args) {
        String folderPath = "D:\\Sem-X\\01_05\\adsorptionIsotherm\\reader"; // Change this to your folder path
        List<FileData> fileDataList = new ArrayList<>();
        //  List<FileDataGraphene> fileDataListGraphene = new ArrayList<>();
        try {
            Files.walk(Paths.get(folderPath))
                    .filter(Files::isRegularFile)
                    .filter(path -> path.toString().endsWith(".txt"))
                    .forEach(path -> {
                        try {
                            List<String> lines = Files.readAllLines(path);
                            String fileName = path.getFileName().toString();
                            String value1 = "";
                            double value2 = 0;
                            String value4 =" ";
                         /*   for (String line : lines) {
                                if (line.startsWith("F://Avagadro")) {
                                    value4 = line.trim();
                                }
                                if (line.startsWith("-")) {
                                    //String[] parts = line.split(" ");
                                    //value1 = (parts[1].trim().replace("+", ""));
                                    value1 = line;
                                }
                                if (line.startsWith("rho ")) {
                                    String[] parts = line.split(" ");
                                    //  String part1 = parts[1].replaceAll("\\d+\\.\\d+", "");
                                    //   String part2 = parts[2].replaceAll("\\(.*\\)", "");
                                    String part1 = parts[1].trim().replace("(", "");
                                    String part2 = parts[2].trim().replace(",", "");
                                 //   String part3 = parts[6].trim().replace(",", "");
                                  //  System.out.println(part1Mod + " " + part2 + " " +part3 +" "+part4);

                                    //Truncate part2 to five decimal places
                                    Double part1Double = Double.parseDouble(part1);
                                    Double part2Double = Double.parseDouble(part2);
                                  //  Double part3Double = Double.parseDouble(part3);
                                    //  String results = removeFirstFourAlphabets(parts[1]);
                                    fileDataList.add(new FileData(fileName, value1, part1Double, part2Double, value4));
                                   // System.exit(2);
                                }
                            }*/
                            for (int i = 0; i < lines.size(); i++) {
                                String line = lines.get(i);

                                if (line.startsWith("F://Avagadro")) {
                                    value4 = line.trim();
                                }

                                if (!line.trim().startsWith("-")) {
                                    value1 = line.trim();  // Store this -3000 line

                                    // Check if there is a line two lines after
                                    if (i  < lines.size()) {
                                        String rhoLine = lines.get(i);
                                        if (rhoLine.startsWith("Num atoms One :")) {
                                            String[] parts = rhoLine.split(" ");
                                            if (parts.length >= 2) {
                                                try {
                                                    String part1 = parts[4].trim().replace("(", "");
                                                    String part2 = parts[5].trim().replace(",", "");
                                                    double part1Double = Double.parseDouble(part1);
                                                    double part2Double = Double.parseDouble(part2);
                                                    //System.out.println(part1Double +" "+ part2Double);
                                                    fileDataList.add(new FileData(fileName, value1, part1Double, part2Double, value4));
                                                } catch (NumberFormatException e) {
                                                    System.err.println("Error parsing rho values in file: " + fileName);
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    });
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("New File");
        for (FileData data : fileDataList) {
            System.out.println(data);
        }
        System.out.println("\n");
    }
}
/*class FileDataGraphene {
    private String fileName;
    private double value1;
    private double value2;

    public FileDataGraphene(String fileName, double value1, double value2) {
        this.fileName = fileName;
        this.value1 = value1;
        this.value2 = value2;
    }

    public String getFileName() {
        return fileName;
    }

    public double getValue1() {
        return value1;
    }

    public double getValue2() {
        return value2;
    }

    @Override
    public String toString() {
        return String.valueOf("value 1" +  value1+"value 2 " + value2 + " value 3 " + value3) ;
    }
}*/

/*public class ReadFiles {
    public static void main(String[] args) {
        String folderPath = "D://SemVI//Data//sorption//0626//sorption//reader"; // Change this to your folder path
        List<FileData> fileDataList = new ArrayList<>();
      //  List<FileDataGraphene> fileDataListGraphene = new ArrayList<>();
        try {
            Files.walk(Paths.get(folderPath))
                    .filter(Files::isRegularFile)
                    .filter(path -> path.toString().endsWith(".txt"))
                    .forEach(path -> {
                        try {
                            List<String> lines = Files.readAllLines(path);
                            String fileName = path.getFileName().toString();
                            double value1 = 0;
                            double value2 = 0;

                            for (String line : lines) {
                                if (line.startsWith("-")) {
                                    value1 = Double.parseDouble(line.trim());
                                }

                                if (line.contains("rho ")) {
                                    String[] parts = line.split(" ");
                                  //  String part1 = parts[1].replaceAll("\\d+\\.\\d+", "");
                                 //   String part2 = parts[2].replaceAll("\\(.*\\)", "");
                                    String part1 = parts[1];
                                    String part2 = parts[2];
                                    //Truncate part2 to five decimal places
                                    Double part1Double = Double.parseDouble(part1);
                                    Double part2Double = Double.parseDouble(part2);
                                    //  String results = removeFirstFourAlphabets(parts[1]);
                                    fileDataList.add(new FileData(fileName, value1, part1Double, part2Double));
                                }
                            }
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    });
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("New File");
        for (FileData data : fileDataList) {
            System.out.println(data);
        }
        System.out.println("\n");
    }
}*/