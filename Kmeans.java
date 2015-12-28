/*******************************************************************************
 * Kmeans.java      Author: Richard Bennett     bennett.cpp@gmail.com
 * ........................................     12/03/2015
 * 
 * input: A file containing selected corrected inputs
 * output: Console
 * 
 * Description: K-Means for Intrusion Detection.
 ******************************************************************************/
import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class Kmeans {
    //Global Variables
    private final static int MAXCLUSTERS = 100;
    private final static int INCREMENT = 20;
    private final static int RUNS = 30;
    //

    //getInput() gets the input from the requested file
    private static String[][] getInput() throws FileNotFoundException, IOException {
        //Variables
        Scanner scan;
        String[][] inputArray; 
        String[] readLine;          
        String inputFile = null;
        int rows;                          
        int cols; 
        File input = null;   
        BufferedReader inputSize;
        BufferedReader inputContents;
        boolean NoFile = true;
        //
        
        //Loop until a valid file is entered
        while (NoFile) {
            scan = new Scanner(System.in);
            System.out.println("Please enter the file name containing the user commands: ");
            inputFile = scan.nextLine();
            input = new File(inputFile);
            if (input.exists()) {
                NoFile = false;
            }
        }
        inputSize = new BufferedReader(new FileReader(input));   
        inputContents = new BufferedReader(new FileReader(input)); 

        //Find the number of cols in the input file
        readLine = inputSize.readLine().split(",");
        cols = readLine.length;
        rows = 1;

        //Find the number of rows in the input file
        while (inputSize.readLine() != null) {
            rows++;
        }
        inputSize.close();

        System.out.println("\nInput file name: " + inputFile);
        System.out.println("* Checking up to groups of " + MAXCLUSTERS + " clusters.");
        System.out.println("* Incrementing " + INCREMENT + " between groups.");
        System.out.println("* For each group refining " + RUNS + " times.\n");

        //Read the input file to the input array
        inputArray = new String[rows][cols];
        for (int i = 0; i < rows; i++) {
            readLine = inputContents.readLine().split(",");
            System.arraycopy(readLine, 0, inputArray[i], 0, cols);
        }
        inputContents.close();
        System.out.println("Input Array initialized:");
        System.out.println("\tInput Array number of rows: " + inputArray.length);
        System.out.println("\tInput Array number of cols: " + inputArray[0].length);
        return inputArray;
    }

    //getData() turns the input into a dataArray
    private static double[][] getData(String[][] inputArray) throws IOException {
        //Variables    
        int rows = inputArray.length;
        int cols = inputArray[0].length - 4; //The number of cols will shrink by four once Strings are removed
        double[][] dataArray = new double[rows][cols];   //holds the data from the input data
        //
        
        //Find all the non-String data in the Input Array
        for (int i = 0; i < rows; i++) {
            dataArray[i][0] = Double.parseDouble(inputArray[i][0]);
            for (int j = 1; j < cols - 4; j++) {
                dataArray[i][j] = Double.parseDouble(inputArray[i][j + 3]);
            }
        }
        System.out.println("Data Array initialized:");
        System.out.println("\tData Array number of rows: " + dataArray.length);
        System.out.println("\tData Array number of cols: " + dataArray[0].length);
        return dataArray;
    }

    //getLabels() makes an array out of all the final strings of each input to use as labels
    private static String[] getLabel(String[][] inputArray) {
        //Variables
        int rows = inputArray.length;
        int cols = inputArray[0].length;
        String[] labelArray = new String[rows];
        //
        
        for (int i = 0; i < rows; i++) {
            labelArray[i] = inputArray[i][cols - 1];
        }
        return labelArray;
    }

    //calcMean() Calculates the mean of each column of the Data Array
    private static double[] calcMean(double[][] dataArray) {
        //Variables
        int rows = dataArray.length;
        int cols = dataArray[0].length;
        double[] mean;      //Holds the mean for each column
        //

        //Add up the sum of each column
        double[] sum = new double[cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sum[j] += dataArray[i][j];
            }
        }
        
        //Divide by the number of rows
        mean = new double[cols];
        for (int i = 0; i < cols; i++) {
            mean[i] = sum[i] / rows;
        }
        return mean;
    }

    //calcStdDev calculates the standard deviation of each column of the Data Array
    private static double[] calcStdDev(double[][] dataArray) {
        //Variables
        int rows = dataArray.length;
        int cols = dataArray[0].length;
        double[] mean = calcMean(dataArray);
        double[] variance;
        double[] stdDev;
        //
        
        //Finding the Variances
        variance = new double[cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                variance[j] += Math.pow((dataArray[i][j] - mean[j]), 2);
            }
        }
        for (int i = 0; i < cols; i++) {
            variance[i] /= rows;
        }

        //Finding the Standard Deviations
        stdDev = new double[cols];
        for (int i = 0; i < cols; i++) {
            stdDev[i] = Math.sqrt(variance[i]);
        }
        return stdDev;
    }

    //calcNormals calculates all normalizations of the Data Array
    private static double[][] calcNormals(double[][] dataArray) {
        //Variables
        int rows = dataArray.length;
        int cols = dataArray[0].length;
        double[] mean = calcMean(dataArray);
        double[] stdDev = calcStdDev(dataArray);
        double[][] normalArray;
        //

        //Creating normalizations array filled with normalizations of the data array 
        normalArray = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (stdDev[j] != 0) {        //If stdDev is zero the normalization is zero to prevent NaNs
                    normalArray[i][j] = (dataArray[i][j] - mean[j]) / stdDev[j]; //(data - mean)/stdDev
                } else {
                    normalArray[i][j] = 0;
                }
            }
        }
        System.out.println("Normalizations generated.");
        return normalArray;
    }

    //selectInitialClusters selects k number of unique rows as initial Cluster Points
    private static double[][] selectInitialClusters(double[][] array, int numClusters) throws IOException {
        //Variables
        boolean selected = false;           //is true when random unique rows have been selected
        int rows = array.length;            //The number of rows in the initial array
        int cols = array[0].length;         //The number of cols in the initial array
        double[] sumArray = new double[rows];     //Holds the sum of each row from the initial array
        int[] selectedClusters;             //The row number of the selected clusters of unique clusters
        Double[] uniqueSums;               //Sums of the unique Rows
        double[][] uniqueArray;             //Array of the unique rows
        double[][] clusterArray = new double[numClusters][cols];    //Holds the final cluster Array
        Set<Double> hashSet = new LinkedHashSet<>();       //hashset used in finding unique rows
        Random random = new Random();                       //For making random numbers
        //
   
        //To find unique rows we will use the sum of each row, and put them in a
        //hash set, then compare that hashset to the initial data to find the rows
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sumArray[i] += array[i][j] * j;
            }
            hashSet.add(sumArray[i]);
        }

        uniqueSums = new Double[hashSet.size()];   //make an array to hold the unique sums
        uniqueSums = hashSet.toArray(uniqueSums);   //Put the hash set in the unique sums array
        uniqueArray = new double[uniqueSums.length][]; //initialize array to hold unique rows

        //Fill the unique Array with unique Rows
        for (int i = 0; i < uniqueSums.length; i++) {
            for (int j = 0; j < rows; j++) {
                if (uniqueSums[i] == sumArray[j]) {
                    uniqueArray[i] = array[i];
                }
            }
        }
        
        //initialize an array to hold the row number of each selected clusters from the unique array.
        //The amount is the amount of clusters desired.
        selectedClusters = new int[numClusters];

        //Fill the selected clusters with a random row ID for a unique row
        for (int i = 0; i < numClusters; i++) {
            selectedClusters[i] = random.nextInt(uniqueArray.length);
        }

        //Make sure none of the random numbers are the same
        while (!selected) {
            int counter = 0;
            for (int i = 0; i < numClusters; i++) {
                for (int j = i + 1; j < numClusters; j++) {
                    if (selectedClusters[i] == selectedClusters[j]) {
                        selectedClusters[j] = random.nextInt(uniqueArray.length);
                        counter++;
                    }
                }
            }
            if (counter == 0) {
                selected = true;
            }
        }

        //Get your selected rows from the unique Array
        for (int i = 0; i < numClusters; i++) {
            clusterArray[i] = uniqueArray[selectedClusters[i]];
        }
        System.out.println("Initial unique clusters chosen.");
        return clusterArray;
    }
    
    //distance() calculates the distance between two points, used in getDistances()
    private static double distance(double[] pointA, double[] pointB) {
        //Variables
        double distance = 0;
        double cols = pointA.length;
        //
        
        for (int i = 0; i < cols; i++) {
            distance += Math.pow((pointA[i] - pointB[i]), 2);
        }
        distance = Math.sqrt(distance);
        return distance;
    }

    //getDistances() finds the distance of every Normalization from each Cluster Point
    private static double[][] getDistances(double[][] normalArray, double[][] clusterArray) {
        //Variables
        int rows = normalArray.length;
        int clusterRows = clusterArray.length;
        double[][] distanceArray = new double[rows][clusterRows];
        //
        
        //Find the distances from each Cluster Point
        for (int i = 0; i < clusterRows; i++) {      
            for (int j = 0; j < rows; j++) {         
                distanceArray[j][i] = distance(clusterArray[i], normalArray[j]);      
            }
        }
        return distanceArray;
    }

    //inCluster() creates an array showing which Cluster each Data Row is in by Distance
    private static int[][] inCluster(double[][] distance) {
        //Variables
        int rows = distance.length;
        int clusterRows = distance[0].length;
        int[][] inCluster = new int[rows][clusterRows];
        //

        //For each Data Row find which Cluster Point is closest
        for (int i = 0; i < rows; i++) {
            double leastDistance = distance[i][0];
            int cluster = 0;
            for (int j = 0; j < clusterRows; j++) {      
                if (distance[i][j] < leastDistance) { 
                    leastDistance = distance[i][j];
                    cluster = j;
                }
            }
            //If it's in that Cluster we'll put a 1 in the corresponding column
            inCluster[i][cluster] = 1;
        }
        return inCluster;
    }

    //kmeansClustering() finds which cluster a row will be in, then forms a new set
    //of Cluster Points based on the averages of each row in the Cluster
    private static double[][] kmeansClustering(double[][] normalArray, double[][] clusterArray) throws IOException {
        //Variables
        int clusterRows = clusterArray.length;      //The number of clusters
        int rows = normalArray.length;                    //The number of rows in the initial array
        int cols = normalArray[0].length;                 //the number of cols in the initial array
        int[] clusterNum = new int[clusterRows];                    //Holds the number of rows in each cluster
        int[][] inCluster;
        double[][] newCluster = new double[clusterRows][cols];      //Holds the new cluster points
        double[][] distance;        //holds the distance of each row from each cluster point
        //

        distance = getDistances(normalArray, clusterArray);
        inCluster = inCluster(distance);

        //Finds how many Data Rows are in each Cluster
        for (int i = 0; i < clusterRows; i++) {
            for (int j = 0; j < rows; j++) {
                if (inCluster[j][i] > 0) {
                    clusterNum[i]++;
                }
            }
        }
        
        //Add all the normalizations for that Cluster
        for (int i = 0; i < clusterRows; i++) {   
            for (int j = 0; j < rows; j++) {             
                if (inCluster[j][i] > 0) {               
                    for (int k = 0; k < cols; k++) {
                        newCluster[i][k] += normalArray[j][k];    
                    }
                }

            }
        }
        
        //Then divide it by the number of entries in that cluster to get the mean
        for (int i = 0; i < clusterRows; i++) {
            for (int j = 0; j < cols; j++) {
                if (clusterNum[i] > 0) {
                    newCluster[i][j] = newCluster[i][j] / clusterNum[i];
                }

            }
        }      
        return newCluster;
    }

    //calcDistanceAverageScores() calculates the average distance score of each Cluster
    private static double calcDistanceAverageScore(double[][] distance) throws IOException {
        //Variables
        int[][] inCluster = inCluster(distance);
        int rows = inCluster.length;
        int cols = inCluster[0].length;
        double distanceAverageScore = 0;
        DecimalFormat df = new DecimalFormat("#.###");
        //
        
        //Calculates the average distance score for each Cluster
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                if (inCluster[j][i] > 0) 
                    distanceAverageScore += distance[j][i];
            }
        }
        distanceAverageScore = distanceAverageScore / rows;
        System.out.println("\tThe average distance score is: " + distanceAverageScore);
        return distanceAverageScore;
    }

    //sortByCluster() sorts lists each data point and which Cluster it is in
    private static int[] sortByCluster(int[][] inCluster) {
        //Variables
        int[] sortByCluster = new int[inCluster.length];
        //
        
        for (int i = 0; i < inCluster.length; i++) {
            for (int j = 0; j < inCluster[0].length; j++) {
                if (inCluster[i][j] == 1)
                    sortByCluster[i] = j;
            }
        }
        return sortByCluster;
    }

    //clusterContents() sorts the Contents of each Cluster by type and
    //prints the results to System.out
    private static void clusterContents(int[] sortedClusters, String[] labelArray, int numClusters) {
        //Variables
        int[] clusterContents;
        int total;
        //

        System.out.println("\n***************************************");
        System.out.println("The best cluster group contained " + numClusters + " clusters.");
        System.out.println("The composition of the best cluster group:");
        System.out.println("(Clusters only displayed if they have contents.)\n");

        //Sort the data rows of each Cluster by types and counts them
        for (int i = 0; i < numClusters; i++) {
            total = 0;
            clusterContents = new int[3];
            for (int j = 0; j < sortedClusters.length; j++) {
                if (sortedClusters[j] == i) {
                    total++;
                    switch (labelArray[j]) {
                        case "normal.":
                            clusterContents[0]++;
                            break;
                        case "snmpgetattack.":
                            clusterContents[1]++;
                            break;
                        case "smurf.":
                            clusterContents[2]++;
                            break;
                        default:
                            break;
                    }
                }
            }

            //If a cluster is not empty print it to System.outs
            if (total != 0) {
                System.out.print("Cluster " + i + 1 + ": ");
                System.out.print("\tnormals: " + clusterContents[0]);
                System.out.print("\t\tsnmpgetAttacks: " + clusterContents[1]);
                System.out.print("\t\tsmurfs: " + clusterContents[2]);
                System.out.println("\t\tTotal Samples: " + total);
            }
        }
    }

    public static void main(String[] args) throws FileNotFoundException, IOException {
        //Variables
        double[][] dataArray;               //the Data
        double[][] normalizationsArray;     //the Normalizations
        double[][] distanceArray;           //the distances from Cluster Points
        double[][] clusterArray;            //the Clusters
        int[][] bestClusterArray = null;    //The Best Cluster
        int[] sortedClusters;               //A list of data rows in the Best Cluster
        int numClusters;                    //Number of Clusters
        double averageDistanceScore;                    //Average Distance Score
        double bestDistanceScore = Double.MAX_VALUE;    //The best Distance Score
        int bestClusterNum = Integer.MAX_VALUE;         //Number of K in the best Cluster
        String[] labelArray;                            //Types of each Data Row
        String[][] inputArray;                          //Input data
        //

        System.out.println("*********************************************\n"
             + "*      K-MEANS FOR INTRUSION DETECTION      *\n"
             + "*          Author: Richard Bennett          *\n"
             + "*********************************************\n");

        inputArray = getInput();
        dataArray = getData(inputArray);
        labelArray = getLabel(inputArray);       
        normalizationsArray = calcNormals(dataArray);
       
        //Starting at twenty to MAXCLUSTERS, incrementing INCREMENT times each time
        //Grab an initial cluster, and refine it RUNS number of times
        for (int k = 20; k <= MAXCLUSTERS; k += INCREMENT) {
            System.out.println("\nUsing " + k + " Clusters*********");
            numClusters = k;
            clusterArray = selectInitialClusters(normalizationsArray, numClusters);
            System.out.print("Refining: ");
            for (int i = 0; i < RUNS; i++) {
                System.out.print((i + 1) + " ");
                clusterArray = kmeansClustering(normalizationsArray, clusterArray);
            }
            
            //Once it's refined get the Distances and Calculate the average
            //distance score
            distanceArray = getDistances(normalizationsArray, clusterArray);
            System.out.println("\nDistance Array created.");
            averageDistanceScore = calcDistanceAverageScore(distanceArray);
            
            //If it's better than the best average score make it the new average
            //score and designate that cluster as the best one
            if (averageDistanceScore < bestDistanceScore) {
                bestDistanceScore = averageDistanceScore;
                bestClusterNum = k;
                bestClusterArray = inCluster(distanceArray);
            }
        }
        System.out.println("\nThe Best number of Clusters was: " + bestClusterNum);
        System.out.println("The average distance score was: " + bestDistanceScore);

        //Once you have the best cluster print out it's contents
        sortedClusters = sortByCluster(bestClusterArray);
        clusterContents(sortedClusters, labelArray, bestClusterNum);
        System.out.println("\n***************************************");
    }
}
