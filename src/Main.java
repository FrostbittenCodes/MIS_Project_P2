
//Quick and dirty java implementation of phase one for CSE 408 project 1


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.PriorityQueue;
import java.util.HashMap;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.EOFException;
import java.util.ArrayList;
import java.util.Comparator;

public class Main 
{

    public static void main(String[] args) 
    {
        //Declarations and Initializations (To defaults if no command line args used to override)
        String fx, fy, fz, fw;
        int[] closestx;
        int[] closesty;
        int[] closestz;
        int[] closestw;
        
        
        //Parse command line arguments
        
        
        //Input files
        fx = "Data/sampledata/X/1.csv";
        fy = "Data/sampledata/Y/1.csv";
        fz = "Data/sampledata/Z/1.csv";
        fw = "Data/sampledata/W/1.csv";

		String in = "quant.csv";
        
        /*
        //Find and list closest sensors (Task 1)
        closestx = SimilarSensor(fx);
        for(int i = 0; i < 20; i++)
            System.out.println("Closest for x," + (i+1) + ": " + closestx[i*3 + 0] + "," + closestx[i*3 + 1] + "," + closestx[i*3 + 2]);
        
        closesty = SimilarSensor(fy);
        for(int i = 0; i < 20; i++)
            System.out.println("Closest for y," + (i+1) + ": " + closesty[i*3 + 0] + "," + closesty[i*3 + 1] + "," + closesty[i*3 + 2]);
        
        closestz = SimilarSensor(fz);
        for(int i = 0; i < 20; i++)
            System.out.println("Closest for z," + (i+1) + ": " + closestz[i*3 + 0] + "," + closestz[i*3 + 1] + "," + closestz[i*3 + 2]);
        
        closestw = SimilarSensor(fw);
        for(int i = 0; i < 20; i++)
            System.out.println("Closest for w," + (i+1) + ": " + closestw[i*3 + 0] + "," + closestw[i*3 + 1] + "," + closestw[i*3 + 2]);
        
        
        
        
        //Predictive coding
        System.out.println(PredictiveCoding(fx, 1, "output.csv", closestx, 0));
        */
        
        //Quantizer
        UniformQuantize(fx, "quant.csv", 16);
        //HashMap<String,String> result = EO1(in,16,"output.csv");
        //double r = DO1(result, 16);
        //HashMap<String,String> result = EO2(in,16,"output.csv");
        //double r = DO2(result,16);
        //Node tree = EO3(in,4,"output.csv");
        //double r = DO3(tree,4);
        HashMap<String, String> result = EO4(in, 16, "output.csv");
        double r = DO4(result, 16);
		
    }
    
    //Given a multivariate time series (in csv form with each line as a different sensor), returns an array
    // containing the most similar sensors based on L1, L2, and correlation distances/similarities
    // e.g. sensor 2 would have L1 in the 3 index, L2 in the 4, and correlation in the 5 index, while sensor
    // 3's values would be in the 6, 7, 8 indices
    //This function satisfies task 1
    public static int[] SimilarSensor(String series)
    {
        int[] array = new int[60];
        double[][] data = ReadData(series);
        
        //Calculate closest sensors and store in array
        for(int k = 0; k < 20; k++)
        {
            double l1 = 1000000;
            double l2 = 1000000;
            double c = 0;
            
            //Calculate the differences/similarities for one sensor
            for(int j = 0; j < 20; j ++)
            {
                //Skip if same sensor
                if(k == j)
                    continue;
                
                double L1 = 0;
                double L2 = 0;
                double avg1 = 0;
                double avg2 = 0;

                //Find L1/L2/avg1/avg2
                for(int i = 0; i < data[k].length; i++)
                {
                    L1 += Math.abs(data[k][i] - data[j][i]);
                    L2 += Math.pow(data[k][i] - data[j][i], 2);
                    avg1 += data[k][i];
                    avg2 += data[j][i];
                }
                L2 = Math.pow(L2, 0.5);
                avg1 = avg1 / data[k].length;
                avg2 = avg2 / data[j].length;
                
                double sum1 = 0;
                double sum2 = 0;
                double sum3 = 0;
                for(int i = 0; i < data[k].length; i++)
                {
                    sum1 += (data[k][i] - avg1)*(data[j][i] - avg2);
                    sum2 += Math.pow(data[k][i] - avg1, 2);
                    sum3 += Math.pow(data[j][i] - avg2, 2);
                }
                double corr = sum1/Math.pow(sum2*sum3,0.5);

                
                if(L1 < l1)
                {
                    array[k*3 + 0] = j; //L1
                    l1 = L1;
                }
                if(L2 < l2)
                {
                    array[k*3 + 1] = j; //L2
                    l2 = L2;
                }
                if(corr > c)
                {
                    array[k*3 + 2] = j; //corr
                    c = corr;
                }
            }
        }
        
        //return
        return array;
    }
    
    //Given a multivariate time series (in csv form with each line as a different sensor), a predictive scheme,
    // an output file location, and depending on the scheme, an array of closest sensors and a similarity measure,
    // calculates and writes the new data set to the output file location. Also returns the total (absolute) prediction error.
    public static double PredictiveCoding(String series, int scheme, String output, int[] closest, int measure)
    // <editor-fold desc="Predictive Coding Functions" defaultstate="collapsed">
    {
        switch(scheme)
        {
            case 1:
                return PC1(series, output);
            case 2:
                return PC2(series, output);
            case 3:
                return PC3(series, output);
            case 4:
                return PC4(series, output);
            case 5:
                return PC5(series, output, closest, measure);
            case 6:
                return PC5(series, output, closest, measure);
            default:
                System.err.println("Invalid scheme");
                return -1;
        }
    }
    
    //Predictive scheme 1 (No change, write the same data to new file)
    public static double PC1(String series, String output)
    {
        double[][] data = ReadData(series);
        if(data == null)
            return -1;
        
        if(WriteData(data, output) == 0)
            return 0;
        else
            return -1;
    }
    
    //Predictive scheme 2
    public static double PC2(String series, String output)
    {
        double[][] data = ReadData(series);
        double[][] data2 = new double[20][data[0].length];
        double error = 0;
        
        if(data == null)
            return -1;
        
        for(int k = 0; k < 20; k++)
        {
            data2[k][0] = data[k][0];
            for(int i = 1; i < data[k].length; i++)
            {
                double predictor = data[k][i-1];
                data2[k][i] = data[k][i] - predictor;
                error += data2[k][i];
            }
        }
        error = Math.abs(error);
        
        if(WriteData(data2, output) == 0)
            return error;
        else
            return -1;
    }
    
    //Predictive scheme 3
    public static double PC3(String series, String output)
    {
        double[][] data = ReadData(series);
        double[][] data2 = new double[20][data[0].length];
        double error = 0;
        
        if(data == null)
            return -1;
        
        for(int k = 0; k < 20; k++)
        {
            data2[k][0] = data[k][0];
            data2[k][1] = data[k][1];
            for(int i = 2; i < data[k].length; i++)
            {
                double predictor = (data[k][i-1] + data[k][i-2])/2.0;
                data2[k][i] = data[k][i] - predictor;
                error += data2[k][i];
            }
        }
        error = Math.abs(error);
        
        if(WriteData(data2, output) == 0)
            return error;
        else
            return -1;
    }
    
    //Predictive scheme 4
    public static double PC4(String series, String output)
    {
        double[][] data = ReadData(series);
        double[][] data2 = new double[20][data[0].length];
        double error = 0;
        
        if(data == null)
            return -1;
        
        for(int k = 0; k < 20; k++)
        {
            data2[k][0] = data[k][0];
            data2[k][1] = data[k][1];
            for(int i = 2; i < data[k].length; i++)
            {
                double predictor = (2.0*data[k][i-1] + data[k][i-2])/3.0;
                data2[k][i] = data[k][i] - predictor;
                error += data2[k][i];
            }
        }
        error = Math.abs(error);
        
        if(WriteData(data2, output) == 0)
            return error;
        else
            return -1;
    }
    
    //Predictive scheme 5
    public static double PC5(String series, String output, int[] closest, int measure)
    {
        double[][] data = ReadData(series);
        double[][] data2 = new double[20][data[0].length];
        double error = 0;
        
        if(data == null)
            return -1;
        
        for(int k = 0; k < 20; k++)
        {
            int similar = closest[k*3 + measure];
            data2[k][0] = data[k][0];
            for(int i = 1; i < data[k].length; i++)
            {
                double predictor = (2.0*data[k][i-1] + data[similar][i-1])/3.0;
                data2[k][i] = data[k][i] - predictor;
                error += data2[k][i];
            }
        }
        error = Math.abs(error);
        
        if(WriteData(data2, output) == 0)
            return error;
        else
            return -1;
    }
    
    //Predictive scheme 6
    public static double PC6(String series, String output, int[] closest, int measure)
    {
        double[][] data = ReadData(series);
        double[][] data2 = new double[20][data[0].length];
        double error = 0;
        
        if(data == null)
            return -1;
        
        for(int k = 0; k < 20; k++)
        {
            int similar = closest[k*3 + measure];
            data2[k][0] = data[k][0];
            for(int i = 1; i < data[k].length; i++)
            {
                double predictor = (data[k][i-1] + data[similar][i-1])/2.0;
                data2[k][i] = data[k][i] - predictor;
                error += data2[k][i];
            }
        }
        error = Math.abs(error);
        
        if(WriteData(data2, output) == 0)
            return error;
        else
            return -1;
    }
    
    // </editor-fold>
    
    //Given a multivariate time series (in csv form with each line as a different sensor), a quantization scheme,
    // an output file location, and a value setting the number of bits to use for levels (levels = 2^r),
    // quantizes the data using the appropriate scheme and writes the new data set to the output file location.
    public static void Quantize(String series, int scheme, String output, int r)
    // <editor-fold desc="Quantizer Functions" >
    {
        switch(scheme)
        {
            case 1:
                UniformQuantize(series, output, r);
                return;
            case 2:
                GaussianQuantize();
                return;
            case 3:
                ReverseGaussianQuantize();
                return;
            default:
                System.err.println("Invalid scheme");
        }
    }
    
    //Uniformly quantizes the data
    public static void UniformQuantize(String series, String output, int r)
    {
        double max =  2;
        double min = -2;
        double levels = Math.pow(2,r);
        double levelsize = ((max-min)/(levels-1));
        double[][] data = ReadData(series);
        
        for(int j = 0; j < 20; j++)
        {
            for(int i = 0; i < data[j].length; i++)
            {
                data[j][i] = ((Math.round((data[j][i] + max)/levelsize))*levelsize)-max;
            }
        }
        WriteData(data, output);
    }
    
    //Quantizes the data according to Gaussian Distributions
    public static void GaussianQuantize()
    {
        
    }
    
    //Quantizes the data according to Reverse Gaussian Distributions
    public static void ReverseGaussianQuantize()
    {
        
    }
    // </editor-fold>
    
    
    public static HashMap<String, String> Encode(String series, String output, int r, int scheme) {
        switch (scheme) {
            case 1:
                return EO1(series, r, output);

            case 2:
//                return E02(series,r,output);

            case 3:
//                return E03(series,output);

            case 4:
                return EO4(series, r, output);

            default:
                break;

        }
        return null;

    }

    public static double Decode(HashMap<String, String> symbolTable, int scheme, int r) {
        switch (scheme) {
            case 1:
                return DO1(symbolTable, r);
            case 2:
            //return DO2(symbolTable,scheme);
            case 3:
            //return DO3(symbolTable,scheme);
            case 4:
                //return DO4(symbolTable,scheme);

            default:
                break;
        }

        return 0;
    }

    public static HashMap<String, String> EO1(String series, int r, String output) {
        double data[][] = ReadData(series);
        HashMap<String, String> symbolTable = new HashMap();
        int symbolCounter = 0;

        try {
            DataOutputStream out = new DataOutputStream(new FileOutputStream("encode"));
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    String sKey = String.valueOf(data[j][i]);
                    String s = symbolTable.get(sKey);
                    if (s == null) {
                        String binStr = Integer.toBinaryString(symbolCounter);
                        while (binStr.length() < r) {
                            // Add leading zeroes back to binary representation                                                   
                            binStr = "0".concat(binStr);
                        }
                        symbolTable.put(sKey, binStr);
                    }
                    symbolCounter++; // keeping track of how many symbols are written
                    // should be 30*20
                }
            }

            int bitCounter = 0;
            int buffer = 0;
            int buffersWritten = 0;
            out.writeInt(symbolCounter); // first 32 bits in message 
            // are total # of symbols
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    String symbol = symbolTable.get(String.valueOf(data[j][i]));
                    if (symbol != null) {
                        char[] c = symbol.toCharArray();
                        for (int k = 0; k < c.length; k++) {
                            if (bitCounter == 32) {
                                out.writeInt(buffer);
                                buffer = 0;
                                bitCounter = 0;
                                buffersWritten++;
                            }
                            buffer <<= 1; // shift to make space for new bit
                            if (c[k] == '0') {
                                buffer |= 0;  // pack a zero
                            } else {
                                buffer |= 1;  // otherwise pack a 1
                            }
                            bitCounter++;
                        }
                    } else {
                        System.out.println("Critical error: symbol not found");
                    }
                }
            }

            // After writing everything, see if the integer is not fully packed
            if (bitCounter != 0) {
                buffer <<= (32 - bitCounter); // move bits all the way to the left of the integer
                out.writeInt(buffer);
                buffersWritten++;
            }

        } catch (FileNotFoundException e) {
            System.out.println("I/O file open failure");
        } catch (IOException e) {
            System.out.println("I/O binary write failure");
        }

        return symbolTable;
    }

    public static double DO1(HashMap<String, String> symbolTable, int r) {
        int symbolCount = 0;
        int symbolsRead = 0;
        int buffer;
        int bitMask = 1 << 31;  //msb bitmask
        int columnCount = 0;
        int bitCounter = 0;

        try {
            DataInputStream dis = new DataInputStream(new FileInputStream("encode"));
            symbolCount = dis.readInt();
            columnCount = symbolCount / 20;
            double[][] data = new double[20][columnCount];
            Integer curInt;
            buffer = dis.readInt();
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < columnCount; i++) {
                    String symbol = "";
                    for (int k = 0; k < r; k++) {
                        if (bitCounter == 32) {
                            buffer = dis.readInt();
                            bitCounter = 0;
                        }

                        if ((buffer & bitMask) == bitMask) // bit is set
                        {
                            symbol = symbol.concat("1");
                        } else // bit is not set
                        {
                            symbol = symbol.concat("0");
                        }
                        buffer <<= 1;
                        bitCounter++;
                    }

                    // Look up symbol in map and get key.
                    // If found, write to array.
                    for (String key : symbolTable.keySet()) {
                        if ((symbolTable.get(key)).equals(symbol)) {
                            System.out.println("Casting to double: " + key);
                            data[j][i] = Double.parseDouble(key);
                            symbolsRead++;
                            break;
                        }
                    }
                }
            }
            WriteData(data, "decode.csv");
        } catch (EOFException e) {
            System.out.println("End of file.");
            System.out.println(symbolsRead + " symbols read.");
        } catch (FileNotFoundException e) {
            System.out.println("I/O file open failure");
        } catch (IOException e) {
            System.out.println("I/O binary read failure");
        }

        return 0;
    }

    public static HashMap<String, String> EO2(String series, int r, String output) {
        double data[][] = ReadData(series);
        HashMap<String, String> symbolTable = new HashMap();
        int symbolCounter = 0;
        int columnCount = data[0].length;

        try {
            DataOutputStream out = new DataOutputStream(new FileOutputStream("encode"));
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    // Fill up hash map.  Cast to string
                    // Then check if value is in map.
                    // Shouldn't use doubles as keys, unsafe.
                    String sKey = String.valueOf(data[j][i]);
                    String s = symbolTable.get(sKey);
                    if (s == null) {
                        String binStr = Integer.toBinaryString(symbolCounter);
                        while (binStr.length() < r) {
                            // Add leading zeroes back to binary representation                                                   
                            binStr = "0".concat(binStr);
                        }
                        symbolTable.put(sKey, binStr);
                        symbolCounter++; // keeping track of how many symbols are written
                    }
                }
            }

            symbolCounter = 0;
            out.writeInt(columnCount);

            // RLE
            int buffer = 0;
            int bitCounter = 0;
            String previousKey = String.valueOf(data[0][0]);
            int runCounter = 0;
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    String curKey = String.valueOf(data[j][i]);
                    if (curKey.equals(previousKey)) {
                        runCounter++;
                        previousKey = curKey;
                    } else {
                        symbolCounter += runCounter;
                        String symbol = symbolTable.get(previousKey);

                        char[] c = symbol.toCharArray();
                        for (int k = 0; k < c.length; k++) {
                            if (bitCounter == 32) {
                                out.writeInt(buffer);
                                buffer = 0;
                                bitCounter = 0;
                            }
                            buffer <<= 1;
                            if (c[k] == '0') {
                                buffer |= 0;
                            } else {
                                buffer |= 1;
                            }
                            bitCounter++;
                        }

			             //now write frequency
                        String freq = Integer.toBinaryString(runCounter);
                        while (freq.length() < 32) {
                            freq = "0".concat(freq);
                        }

                        char[] z = freq.toCharArray();
                        for (int l = 0; l < z.length; l++) {
                            if (bitCounter == 32) {
                                out.writeInt(buffer);
                                buffer = 0;
                                bitCounter = 0;
                            }
                            buffer <<= 1;
                            if (z[l] == '0') {
                                buffer |= 0;
                            } else {
                                buffer |= 1;
                            }
                            bitCounter++;
                        }
                        runCounter = 1;
                        previousKey = curKey;
                    }

                }
            }

            if (runCounter != 0) {
                symbolCounter += runCounter;
                String symbol = symbolTable.get(previousKey);

                char[] c = symbol.toCharArray();
                for (int k = 0; k < c.length; k++) {
                    if (bitCounter == 32) {
                        out.writeInt(buffer);
                        buffer = 0;
                        bitCounter = 0;
                    }
                    buffer <<= 1;
                    if (c[k] == '0') {
                        buffer |= 0;
                    } else {
                        buffer |= 1;
                    }
                    bitCounter++;
                }

				//now write frequency
                String freq = Integer.toBinaryString(runCounter);
                while (freq.length() < 32) {
                    freq = "0".concat(freq);
                }

                char z[] = freq.toCharArray();
                for (int l = 0; l < z.length; l++) {
                    if (bitCounter == 32) {
                        out.writeInt(buffer);
                        buffer = 0;
                        bitCounter = 0;
                    }
                    buffer <<= 1;
                    if (z[l] == '0') {
                        buffer |= 0;
                    } else {
                        buffer |= 1;
                    }
                    bitCounter++;
                }
            }
            if (bitCounter != 0) {
                buffer <<= (32 - bitCounter); // move bits all the way to the left of the integer
                out.writeInt(buffer);
            }
            System.out.println("Number of symbols read: " + symbolCounter);

        } catch (FileNotFoundException e) {
            System.out.println("I/O file open failure");
        } catch (IOException e) {
            System.out.println("I/O binary write failure");
        }

        return symbolTable;
    }

    public static double DO2(HashMap<String, String> symbolTable, int r) {
        int buffer;
        int bitMask = 1 << 31;
        int columnCount = 0;
        int bitCounter = 0;
        double[][] data = null;
        ArrayList<String> resultList = new ArrayList<String>();

        try {
            DataInputStream dis = new DataInputStream(new FileInputStream("encode"));
            columnCount = dis.readInt();

            data = new double[20][columnCount];
            buffer = dis.readInt();
            while (true) {
                // first read the code
                String symbol = "";
                for (int k = 0; k < r; k++) {
                    if (bitCounter == 32) {
                        buffer = dis.readInt();
                        bitCounter = 0;
                    }

                    if ((buffer & bitMask) == bitMask) // bit is set
                    {
                        symbol = symbol.concat("1");
                    } else // bit is not set
                    {
                        symbol = symbol.concat("0");
                    }
                    buffer <<= 1;
                    bitCounter++;
                }

                //	lookup key
                for (String key : symbolTable.keySet()) {
                    if ((symbolTable.get(key)).equals(symbol)) {
                        //found key, now get frequency
                        String freqString = "";
                        for (int l = 0; l < 32; l++) {
                            if (bitCounter == 32) {
                                buffer = dis.readInt();
                                bitCounter = 0;
                            }
                            if ((buffer & bitMask) == bitMask) {
                                freqString = freqString.concat("1");
                            } else {
                                freqString = freqString.concat("0");
                            }
                            buffer <<= 1;
                            bitCounter++;
                        }
                        int finalFreq = Integer.parseInt(freqString, 2); //parse binary to int
                        for (int x = 0; x < finalFreq; x++) {
                            resultList.add(key);
                        }
                        break;
                    }
                }
            }
        } catch (EOFException e) {
            System.out.println("Size of list: " + resultList.size());
            int curInt = 0;
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < columnCount; i++) {
                    data[j][i] = Double.parseDouble(resultList.get(curInt));
                    curInt++;
                }
            }
            WriteData(data, "decode.csv");

        } catch (FileNotFoundException e) {
            System.out.println("I/O file open failure");
        } catch (IOException e) {
            System.out.println("I/O binary write failure");
        }
        return 0;
    }

    public static Node EO3(String series, int r, String output) {
        double data[][] = ReadData(series);
        HashMap<String, Integer> freqTable = new HashMap();
        HashMap<String, String> symbolTable = new HashMap();
        int columnCount = data[0].length;
        PriorityQueue<Node> q = new PriorityQueue<Node>(10, nodeComparator);
        Node tree = null;

        try {
            DataOutputStream out = new DataOutputStream(new FileOutputStream("encode"));
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    String sKey = String.valueOf(data[j][i]);
                    Integer s = freqTable.get(sKey);
                    if (s == null) {
                        freqTable.put(sKey, 0);
                    } else {
                        freqTable.put(sKey, freqTable.get(sKey) + 1);
                    }
                }
            }

            for (String key : freqTable.keySet()) {
                Node newNode = new Node();
                newNode.symbol = key;
                newNode.value = freqTable.get(key);
                q.add(newNode);
            }

            Node root = buildTree(q);

            PreOrder(root, "");
            FillTable(root, symbolTable);

            tree = root;

            //Now we have the symbols... write to file
            int bitCounter = 0;
            int buffer = 0;
            out.writeInt(columnCount);
            int buffersWritten = 0;
            Boolean first = true;
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    String symbol = symbolTable.get(String.valueOf(data[j][i]));
                    if (symbol != null) {
                        char[] c = symbol.toCharArray();
                        for (int k = 0; k < c.length; k++) {
                            if (bitCounter == 32) {
                                out.writeInt(buffer);
                                buffer = 0;
                                bitCounter = 0;
                                buffersWritten++;
                            }
                            buffer <<= 1;
                            if (c[k] == '0') {
                                buffer |= 0;
                            } else {
                                buffer |= 1;
                            }
                            bitCounter++;
                        }
                    } else {

                        System.out.println("Critical error symbol not found for " + data[j][i]);
                    }
                }
            }

            if (bitCounter != 0) {
                buffer <<= (32 - bitCounter); // move bits all the way to the left of the integer
                out.writeInt(buffer);
                buffersWritten++;
            }
            System.out.println("Buffers written: " + buffersWritten);

        } catch (FileNotFoundException f) {
            System.out.println("Error");
        } catch (IOException e) {
            System.out.println("I/O Binary write failure");
        }

        return tree;
    }

    public static double DO3(Node head, int r) {
        int buffer;
        int bitMask = 1 << 31;
        int columnCount = 0;
        int bitCounter = 0;
        double[][] data = null;
        int printCounter = 0;
        Boolean first = true;

        try {
            DataInputStream dis = new DataInputStream(new FileInputStream("encode"));
            columnCount = dis.readInt();
            data = new double[20][columnCount];
            buffer = dis.readInt();

            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    Node root = head;
                    Node prev = null;
                    while (root.right != null || root.left != null) {
                        String finalSymbol = "";
                        String finalKey = "";
                        if (bitCounter == 32) {
                            buffer = dis.readInt();
                            bitCounter = 0;
                        }

                        if ((buffer & bitMask) == bitMask) {
                            prev = root;
                            root = root.right;
                            buffer <<= 1;
                            bitCounter++;

                        } else {
                            prev = root;
                            root = root.left;
                            buffer <<= 1;
                            bitCounter++;
                        }
                    }
                    first = false;
                    data[j][i] = Double.parseDouble(root.symbol);
                    // System.out.println("Data is: " + data[j][i]);
                }
            }

        } catch (EOFException e) {
            System.out.println("End of file.");
        } catch (FileNotFoundException e) {
            System.out.println("I/O file open failure");
        } catch (IOException e) {
            System.out.println("I/O binary read failure");
        }
        WriteData(data, "decode.csv");

        return 0;

    }

    public static HashMap<String, String> EO4(String series, int r, String output) {
        double data[][] = ReadData(series);
        HashMap<String, String> symbolTable = new HashMap();
        int symbolCounter = 0;
        int codeCounter = 0;
        int buffersWritten = 0;
        int columnCount = data[0].length;
        ArrayList<String> inputList = new ArrayList<String>();

        try {
            // Fill up table with primary symbols first
            DataOutputStream out = new DataOutputStream(new FileOutputStream("encode"));
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    String sKey = String.valueOf(data[j][i]);
                    String s = symbolTable.get(sKey);
                    if (s == null) {
                        String codeStr = Integer.toBinaryString(codeCounter);
                        while (codeStr.length() < r) {
                            codeStr = "0".concat(codeStr);
                        }
                        symbolTable.put(sKey, codeStr);
                        codeCounter++;
                    }
                    inputList.add(sKey);
                }
            }

            // LZW compression
            int buffer = 0;
            int bitCounter = 0;

            // First write num of columns to file
            out.writeInt(columnCount);

            String s = inputList.get(0);
            String c = null;

            for (int i = 1; i < inputList.size(); i++) {
                c = inputList.get(i);
                if ((symbolTable.get(s.concat(",").concat(c))) != null) {
                    s = s.concat(",").concat(c);
                } else {
                    // Output code for s
                    String code = symbolTable.get(s);
                    char[] ch = code.toCharArray();
                    for (int k = 0; k < ch.length; k++) {
                        if (bitCounter == 32) {
                            out.writeInt(buffer);
                            buffer = 0;
                            bitCounter = 0;
                        }

                        buffer <<= 1;
                        if (ch[k] == '0') {
                            buffer |= 0;
                        } else {
                            buffer |= 1;
                        }
                        bitCounter++;
                    }

                    // Add s+c to dict with new code
                    String codeStr = Integer.toBinaryString(codeCounter);
                    while (codeStr.length() < r) {
                        codeStr = "0".concat(codeStr);
                    }
                    symbolTable.put(s.concat(",").concat(c), codeStr);
                    codeCounter++;

                    // s = c
                    s = c;
                }
            }

            // Output code for s;
            String code = symbolTable.get(s);
            char[] ch = code.toCharArray();
            for (int k = 0; k < ch.length; k++) {
                if (bitCounter == 32) {
                    out.writeInt(buffer);
                    buffer = 0;
                    bitCounter = 0;
                }

                buffer <<= 1;
                if (ch[k] == '0') {
                    buffer |= 0;
                } else {
                    buffer |= 1;
                }
                bitCounter++;
            }

            // Clear buffer if there's bits leftover
            if (bitCounter != 0) {
                buffer <<= (32 - bitCounter); // move bits all the way to the left of the integer
                out.writeInt(buffer);
            }

        } catch (FileNotFoundException e) {
            System.out.println("I/O file open failure");
        } catch (IOException e) {
            System.out.println("I/O binary write failure");
        }
        return symbolTable;
    }

    public static double DO4(HashMap<String, String> symbolTable, int r) {
        int buffer;
        int bitMask = 1 << 31;
        int columnCount = 0;
        int bitCounter = 0;
        String finalResult = "";

        try {
            DataInputStream dis = new DataInputStream(new FileInputStream("encode"));
            columnCount = dis.readInt();
            buffer = dis.readInt();
            while (true) {
                String symbol = "";

                for (int i = 0; i < r; i++) {
                    if (bitCounter == 32) {
                        buffer = dis.readInt();
                        bitCounter = 0;
                    }

                    if ((buffer & bitMask) == bitMask) {
                        symbol = symbol.concat("1");
                    } else {
                        symbol = symbol.concat("0");
                    }

                    buffer <<= 1;
                    bitCounter++;
                }
                for (String key : symbolTable.keySet()) {

                    if ((symbolTable.get(key)).equals(symbol)) {
                        finalResult = finalResult.concat(",").concat(key);
                        break;
                    }
                }
            }

        } catch (EOFException e) {
            System.out.println("End of file");
            String[] valueList = finalResult.split(",");

            double[][] data = new double[20][columnCount];
            int current = 1; // Start at 2nd value - first value is always null.

            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < columnCount; i++) {
                    if (current < valueList.length) {
                        data[j][i] = Double.parseDouble(valueList[current]);
                        current++;
                    }

                }
            }
            WriteData(data, "decode.csv");

        } catch (FileNotFoundException e) {
            System.out.println("I/O file open failure");
        } catch (IOException e) {
            System.out.println("I/O binary read failure");
        }

        return 0;
    }
    //Everything below this line is a supporting function (e.g. read/write arrays from/to files)
    //------------------------------------------------------------------------------------------
    // <editor-fold desc="Supporting Functions" defaultstate="collapsed">
    
    //Reads data in from a csv and returns an 2d array of doubles with 20 rows
    public static double[][] ReadData(String series)
    {
        double[][] data = new double[20][];
        BufferedReader input;
        
        try
        {
            input = new BufferedReader(new FileReader(series));
        }
        //Error occurred, return null
        catch(FileNotFoundException e)
        {
            System.err.println(e);
            return null;
        }
        
        //Read and parse the data and to populate the data array
        try
        {
            //Executes 20 times (once for each sensor)
            for(int j = 0; j < 20; j++)
            {
                String line;
                //Read and parse data for the sensors
                line = input.readLine();
                String temp[] = line.split(",");
                data[j] = new double[temp.length];
                for(int i = 0; i < temp.length; i++)
                {
                    data[j][i] = Double.parseDouble(temp[i]);
                }
            }
        }
        //Error occurred, return null
        catch(IOException e)
        {
            System.err.println(e);
            return null;
        }
        
        try
        {
            input.close();
        }
        catch(IOException e)
        {
            System.err.println(e);
            return null;
        }
        
        
        return data;
    }
    
    //Writes a 2d array to a csv file, returns 0 on success, -1 on failure
    public static int WriteData(double[][] data, String filepath)
    {
        BufferedWriter output;
        String line;
        
        try
        {
            output = new BufferedWriter(new FileWriter(filepath));
        }
        //Error occurred, return -1
        catch(IOException e)
        {
            System.err.println(e);
            return -1;
        }
        
        try
        {
            for(int j = 0; j < 20; j++)
            {
                line = "";
                for(int i = 0; i < data[j].length; i++)
                {
                    line += data[j][i];
                    line += ',';
                }
                line += '\n';
                output.write(line);
            }
        }
        catch(IOException e)
        {
            System.err.println(e);
            return -1;
        }
        
        try
        {
            output.flush();
            output.close();
        }
        catch(IOException e)
        {
            System.err.println(e);
            return -1;
        }
        
        return 0;
    }

    public static Node buildTree(PriorityQueue<Node> q) {
        int parentSymbol = 0;
        Node parent = null;
        Node a = null;
        Node b = null;

        if (q.size() < 2) {
            return null;
        } else {
            while (q.size() > 1) {
                a = q.poll();
                b = q.poll();
                parent = new Node();
                parent.value = a.value + b.value;
                parent.symbol = "P".concat(Integer.toString(parentSymbol));
                a.parent = parent;
                b.parent = parent;
                parent.right = b;
                parent.left = a;
                parentSymbol++;
                q.add(parent);
            }
            return q.poll();
        }
    }

    public static void PreOrder(Node n, String id) {
        if (n == null) {
            return;
        }
        n.code = id;
        PreOrder(n.left, id.concat("0"));
        PreOrder(n.right, id.concat("1"));
    }

    public static void FillTable(Node n, HashMap<String, String> symbolTable) {
        if (n == null) {
            return;
        }
        if (n.code != null) {
            symbolTable.put(n.symbol, n.code);
        }
        FillTable(n.left, symbolTable);
        FillTable(n.right, symbolTable);
    }

    public static Comparator<Node> nodeComparator = new Comparator<Node>() {
        @Override
        public int compare(Node n1, Node n2) {
            return (int) (n1.value - n2.value);
        }
    };


    // </editor-fold>
}
