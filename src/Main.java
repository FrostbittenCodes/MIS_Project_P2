
//Quick and dirty java implementation of phase one for CSE 408 project 1


import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import javax.imageio.ImageIO;

public class Main 
{

    public static void main(String[] args) 
    {
        //Declarations and Initializations (To defaults if no command line args used to override)
        String fx, fy, fz, fw, output_root, root;
        String fxn, fyn, fzn, fwn;
        int encode_scheme = 1;
        int quant_scheme = 1;
        int pc_scheme = 0;
        int closest_scheme = 0;
        int bits = 8;
        double ex = 0, ey = 0, ez = 0, ew = 0;
        int[] closestx;
        int[] closesty;
        int[] closestz;
        int[] closestw;
        
        
        //Parse command line arguments
        for(int i = 0; i < args.length-1; i++)
        {
            switch(args[i])
            {
                //Help flag
                case"-h":
                    System.out.println("Usage: java -jar Main.jar [options] <infile root>\nSee readme for options.");
                    return;
                //Quantizer scheme file flag
                case "-b":
                    bits = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                //Quantizer scheme file flag
                case "-q":
                    quant_scheme = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                //Encoder scheme file flag
                case "-e":
                    encode_scheme = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                //Closest scheme file flag
                case "-c":
                    closest_scheme = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                //Prediction scheme file flag
                case "-p":
                    pc_scheme = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                default:
                    System.err.println("Unrecognized flag");
                    return;
            }
        }
        //Path to root
        root = args[args.length-1];
        
        
        //Input files
        fx = root + "/X/1.csv";
        fy = root + "/Y/1.csv";
        fz = root + "/Z/1.csv";
        fw = root + "/W/1.csv";
        
        fxn = fx;
        fyn = fy;
        fzn = fz;
        fwn = fw;
        
        //set output root
        output_root = root + "/output/";        
        
        
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
        
        
        //Predictive coding (Task 2
        if(pc_scheme == 0)
            System.out.println("Using no pc");
        else if(pc_scheme > 0 && pc_scheme < 5)
        {
            fxn = output_root + "xpc.csv";
            fyn = output_root + "ypc.csv";
            fzn = output_root + "zpc.csv";
            fwn = output_root + "wpc.csv";
            ex = PredictiveCoding(fx, pc_scheme, fxn, closestx, 0);
            ey = PredictiveCoding(fy, pc_scheme, fyn, closesty, 0);
            ez = PredictiveCoding(fz, pc_scheme, fzn, closestz, 0);
            ew = PredictiveCoding(fw, pc_scheme, fwn, closestw, 0);
        }
        else if(pc_scheme >= 5 && pc_scheme < 7)
        {
            fxn = output_root + "xpc.csv";
            fyn = output_root + "ypc.csv";
            fzn = output_root + "zpc.csv";
            fwn = output_root + "wpc.csv";
            if(closest_scheme >= 0 && closest_scheme < 3)
            {
                ex = PredictiveCoding(fx, pc_scheme, fxn, closestx, closest_scheme);
                ey = PredictiveCoding(fy, pc_scheme, fyn, closesty, closest_scheme);
                ez = PredictiveCoding(fz, pc_scheme, fzn, closestz, closest_scheme);
                ew = PredictiveCoding(fw, pc_scheme, fwn, closestw, closest_scheme);
            }
            else
            {
                System.err.println("Unrecognized Similarity Scheme, using L1 instead");
                ex = PredictiveCoding(fx, pc_scheme, fxn, closestx, 0);
                ey = PredictiveCoding(fy, pc_scheme, fyn, closesty, 0);
                ez = PredictiveCoding(fz, pc_scheme, fzn, closestz, 0);
                ew = PredictiveCoding(fw, pc_scheme, fwn, closestw, 0);
            }
        }
        else
            System.err.println("Unrecognized Predictive Scheme, using no pc instead");
        
        System.out.println("Absolute error resulting from Predictive Coding:");
        System.out.println("x: " + ex);
        System.out.println("y: " + ey);
        System.out.println("z: " + ez);
        System.out.println("w: " + ew);
        
        //Point to new location if necessary
        fx = fxn;
        fy = fyn;
        fz = fzn;
        fw = fwn;
        
        
        //Quantizer (Task 3)
        if(quant_scheme == 0)
            System.out.println("Using no quantization");
        else if(quant_scheme > 0 && quant_scheme < 4)
        {
            fxn = output_root + "xquant.csv";
            fyn = output_root + "yquant.csv";
            fzn = output_root + "zquant.csv";
            fwn = output_root + "wquant.csv";
            Quantize(fx, quant_scheme, fxn, bits);
            Quantize(fy, quant_scheme, fyn, bits);
            Quantize(fz, quant_scheme, fzn, bits);
            Quantize(fw, quant_scheme, fwn, bits);
        }
        else
            System.err.println("Unrecognized Quantization Scheme, using no quantization instead");
        
        fx = fxn;
        fy = fyn;
        fz = fzn;
        fw = fwn;
        
        //Encode (Task 4)
        
        //Decode (Task 5)
        
        //Visualizer
        VisualizeNoise(root + "/X/1.csv", fx, root + "xdiff.bmp");
        VisualizeNoise(root + "/Y/1.csv", fy, root + "ydiff.bmp");
        VisualizeNoise(root + "/Z/1.csv", fz, root + "zdiff.bmp");
        VisualizeNoise(root + "/W/1.csv", fw, root + "wdiff.bmp");
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
                GaussianQuantize(series, output, r);
                return;
            case 3:
                ReverseGaussianQuantize(series, output, r);
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
        double[][] lookup = new double[2][(int)levels];
        
        //populate lookup table entries
        lookup[0][0] = Double.NEGATIVE_INFINITY;
        lookup[0][1] = min+levelsize/2;
        for(int i = 0; i < levels; i++)
        {
            if(i > 1)
                lookup[0][i] = lookup[0][i-1]+levelsize;
            lookup[1][i] = min+i*levelsize;
        }
        
        //Quantize using lookup table
        for(int k = 0; k < 20; k++)
        {
            for(int j = 0; j < data[k].length; j++)
            {
                for(int i = 0; i < levels-1; i++)
                {
                    if(data[k][j] >= lookup[0][i+1]){}
                    else
                    {
                        data[k][j] = lookup[1][i];
                        break;
                    }
                }
            }
        }
        WriteData(data, output);
    }
    
    //Quantizes the data according to Gaussian Distributions
    public static void GaussianQuantize(String series, String output, int r)
    {
        double max =  2;
        double min = -2;
        double levels = Math.pow(2,r);
        double levelsize = ((max-min)/(levels-1));
        double[][] data = ReadData(series);
        double[][] lookup = new double[2][(int)levels];
        double[] tempvals = new double[(int)levels];
        
        //Equal bands
        for(int i = 0; i < levels; i++)
        {
            tempvals[i] = min+i*levelsize;
        }
        
        //Gaussian bands (mirrored it to make it symmetrical
        // if not mirrored, fp errors do strange things to it
        double sum = min;
        
        lookup[1][0] = -2;
        lookup[0][0] = Double.NEGATIVE_INFINITY;
        for(int i = 1; i < levels/2; i++)
        {
            double temp = cdf(tempvals[i]) - cdf(tempvals[i-1]);
            sum+=temp*(max-min);
            lookup[1][i] = sum;
            lookup[1][(int)levels-i] = -lookup[1][i-1];
        }
        lookup[1][(int)levels/2] = -lookup[1][(int)levels/2-1];
        
        
        //Populate the rest of the lookup table
        for(int i = 1; i < levels; i++)
        {
            if(i > 1)
                lookup[0][i] = (lookup[1][i-1] + lookup[1][i])/2;
            else
                lookup[0][i] = (min + lookup[1][i])/2;
        }
        
        
        //Quantize using lookup table
        for(int k = 0; k < 20; k++)
        {
            for(int j = 0; j < data[k].length; j++)
            {
                for(int i = 0; i < levels-1; i++)
                {
                    if(data[k][j] >= lookup[0][i+1]){}
                    else
                    {
                        data[k][j] = lookup[1][i];
                        break;
                    }
                }
            }
        }
        WriteData(data, output);
    }
    
    //Quantizes the data according to Reverse Gaussian Distributions
    public static void ReverseGaussianQuantize(String series, String output, int r)
    {
        double max =  2;
        double min = -2;
        double levels = Math.pow(2,r);
        double levelsize = ((max-min)/(levels-1));
        double[][] data = ReadData(series);
        double[][] lookup = new double[2][(int)levels];
        double[] tempvals = new double[(int)levels];
        
        //Equal bands
        for(int i = 0; i < levels; i++)
        {
            tempvals[i] = min+i*levelsize;
        }
        
        //Gaussian bands (mirrored it to make it symmetrical
        // if not mirrored, fp errors make the positives slightly different
        double sum = 0;
        
        lookup[1][0] = -2;
        lookup[0][0] = Double.NEGATIVE_INFINITY;
        for(int i = 1; i < levels/2; i++)
        {
            double temp = cdf(tempvals[i]) - cdf(tempvals[i-1]);
            sum+=temp*(min-max);
            lookup[1][i] = sum;
        }

        //Reverse the generated stretch
        for(int i = 1; i < levels/4; i++)
        {
            double temp = lookup[1][i];
            lookup[1][i] = lookup[1][(int)levels/2-i];
            lookup[1][(int)levels/2-i] = temp;
        }
        
        //Mirror
        for(int i = 1; i < levels/2; i++)
        {
            lookup[1][(int)levels-i] = -lookup[1][i-1];
        }
        lookup[1][(int)levels/2] = -lookup[1][(int)levels/2-1];
        
        //Populate the rest of the lookup table
        for(int i = 1; i < levels; i++)
        {
            if(i > 1)
                lookup[0][i] = (lookup[1][i-1] + lookup[1][i])/2;
            else
                lookup[0][i] = (min + lookup[1][i])/2;
        }
        
        //Quantize using lookup table
        for(int k = 0; k < 20; k++)
        {
            for(int j = 0; j < data[k].length; j++)
            {
                for(int i = 0; i < levels-1; i++)
                {
                    if(data[k][j] >= lookup[0][i+1]){}
                    else
                    {
                        data[k][j] = lookup[1][i];
                        break;
                    }
                }
            }
        }
        WriteData(data, output);
    }
    // </editor-fold>
    
    //Given two multivariate time series (in csv form with each line as a different sensor), outputs an 
    // image showing the difference between the two input files.
    // Currently doesn't calculate noise, but outputs the appropriate visual stuff
    public static void VisualizeNoise(String f1, String f2, String out)
    // <editor-fold>
    {
        BufferedReader input1;
        BufferedReader input2;
        String outputfile = out;
        double[][] data1 = new double[20][];
        double[][] data2 = new double[20][];
        double[][] diff = new double[20][];
        int B = 1; //Beta value
        int maxlength;
        int[] imagedata;
        int width;
        int padding = 20;
        int thickness = 10;
        int height = thickness*20 + 22*padding;
        int color1 = 0x0000ff; //lower color
        int color2 = 0xff0000; //upper color
        BufferedImage output;
        int background = 0xffffff;
        
        
        //Open the two file readers
        try
        {
            input1 = new BufferedReader(new FileReader(f1));
            input2 = new BufferedReader(new FileReader(f2));
        }
        catch(FileNotFoundException e)
        {
            System.err.println(e);
            return;
        }
        
        //Read and parse the data for files
        try
        {
            String line;
            //Executes 20 times (for the given data sets)
            for(int j = 0; j < 20; j++)
            {
                //Read and parse data for file 1
                line = input1.readLine();
                String temp[] = line.split(",");
                data1[j] = new double[temp.length];
                for(int i = 0; i < temp.length; i++)
                {
                    data1[j][i] = Double.parseDouble(temp[i]);
                }
                //Read and parse data for file 2
                line = input2.readLine();
                temp = line.split(",");
                data2[j] = new double[temp.length];
                for(int i = 0; i < temp.length; i++)
                {
                    data2[j][i] = Double.parseDouble(temp[i]);
                }
            }
        }
        catch(IOException e)
        {
            System.err.println(e);
        }
        
        //Figure out which data set is longer, and store maxlength
        if(data1[0].length > data1[0].length)
            maxlength = data1[0].length;
        else
            maxlength = data2[0].length;
        
        //Calculate the difference and perform scaling (outer loop runs 20 times, inner runs variable (maxlength) amount)
        for(int j = 0; j < 20; j++)
        {
            double max = 0;
            diff[j] = new double[maxlength];
            for(int i = 0; i < maxlength; i++)
            {
                if((data1[j].length > i) && (data2[j].length > i))
                {
                    diff[j][i] = Math.abs(data1[j][i]-data2[j][i]);
                }
                else if(data1[j].length > i)
                {
                    diff[j][i] = Math.abs(data1[j][i]);
                }
                else
                {
                    diff[j][i] = Math.abs(data2[j][i]);
                }
                if(diff[j][i] > max)
                    max = diff[j][i];
            }
            //scale by max in sensor set
            for(int i = 0; i < maxlength; i++)
            {
                diff[j][i] = diff[j][i]/max;
                diff[j][i] = Math.pow(diff[j][i], B);
            }
        }
        
        int temp[] = new int[maxlength*20];
        for(int j = 0; j < 20; j++)
        {
            for(int i = 0; i < diff[j].length; i++)
            {
                temp[(maxlength*j)+i] = getRGBColor(color1, color2, diff[j][i]);
            }
        }
        //Set image width and instantiate imagedata/image
        width = maxlength*thickness + 2*padding;
        imagedata = new int[width*height];
        output = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        
        //Draw background then populate imagedata
        for(int i = 0; i < imagedata.length; i++)
            imagedata[i] = background;
        for(int j = 0; j < 20; j++)
        {
            int r,q;
            for(int i = 0; i < diff[j].length; i++)
            {
                //Sets integers in appropriate width x width area
                // to the specified color rgb value calculated from above
                for(r = 0; r < thickness; r++)
                    for(q = 0; q < thickness; q++)
                        imagedata[(padding + j*(thickness+padding)+r)*width + (padding+i*thickness + q)] = temp[(maxlength*j)+i];
            }
        }
        
        //Create image from data and write to file
        output.setRGB(0, 0, width, height, imagedata, 0, width);
        try
        {
            ImageIO.write(output, "bmp", new File(outputfile));
            System.out.println("Image successfully generated");
        }
        catch(IOException e)
        {
            System.err.println(e);
        }
    }
    // </editor-fold>
    
    
    public static HashMap<String, String> Encode(String series, String output, int r, int scheme) 
    {
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

    public static double Decode(HashMap<String, String> symbolTable, int scheme, int r) 
    {
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

    public static HashMap<String, String> EO1(String series, int r, String output) 
    {
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

            System.out.println("Symbols written: " + symbolCounter);
            System.out.println("Buffers written: " + buffersWritten);
        } catch (FileNotFoundException e) {
            System.out.println("I/O file open failure");
        } catch (IOException e) {
            System.out.println("I/O binary write failure");
        }

        return symbolTable;
    }

    public static double DO1(HashMap<String, String> symbolTable, int r) 
    {
        int symbolCount;
        int symbolsRead = 0;
        int buffer;
        int bitMask = 1 << 31;  //msb bitmask
        int columnCount;
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
                            data[j][i] = Double.parseDouble(key);
                            symbolsRead++;
                            break;
                        }
                    }
                }
            }
            System.out.println("Symbols read: " + symbolsRead);
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

    public static HashMap<String, String> EO2(String series, int r, String output) 
    {
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
                        System.out.println("Run of " + runCounter + " symbols");
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
                System.out.println("Last run is " + runCounter + " symbols");
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

    public static double DO2(HashMap<String, String> symbolTable, int r) 
    {
        int buffer;
        int bitMask = 1 << 31;
        int columnCount = 0;
        int bitCounter = 0;
        double[][] data = null;
        ArrayList<String> resultList = new ArrayList<>();

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
//                        data[j][i] = Double.parseDouble(key);
//                        symbolsRead++; 
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

    public static HashMap<String, String> EO3(String series, int r, String output)
    {
        double data[][] = ReadData(series);
        HashMap<String, Integer> freqTable = new HashMap();
        HashMap<String, String> symbolTable = new HashMap();
        int columnCount = data[0].length;
        PriorityQueue<Node> q = new PriorityQueue<>(10,nodeComparator);
		

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

//			PriorityQueue<Node> temp = new PriorityQueue<Node>(q);
/*			DEBUG PRINTING:
             while (!temp.isEmpty()){
             Node cur = temp.remove();
             System.out.println("Symbol: " + cur.symbol + " Value: " + cur.value);
             }*/
            Node root = buildTree(q);

			// fill up symbol table
            symbolTable.put(root.left.symbol, "0");
            symbolTable.put(root.right.symbol, "1");
            root = root.right;
            while (root.right != null) {
                symbolTable.put(root.left.symbol, symbolTable.get(root.symbol).concat("0"));
                symbolTable.put(root.right.symbol, symbolTable.get(root.symbol).concat("1"));
                root = root.right;
            }

            for (String key : symbolTable.keySet()) {
                System.out.println("Key: " + key + " Value: " + symbolTable.get(key));

            }

        } catch (FileNotFoundException f) {
            System.out.println("Error");
        }
        /*catch(IOException e){
         System.out.println("I/O Binary write failure");
         }*/

        return symbolTable;
    }

    public static HashMap<String, String> EO4(String series, int r, String output) 
    {
        double data[][] = ReadData(series);
        HashMap<String, String> symbolTable = new HashMap();
        int symbolCounter = 0;
        int codeCounter = 0;
        int buffersWritten = 0;
        int columnCount = data[0].length;

        try {
            // Fill up table with primary symbols first
            DataOutputStream out = new DataOutputStream(new FileOutputStream("encode"));
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    String sKey = String.valueOf(data[j][i]);
                    String s = symbolTable.get(sKey);
                    if (s == null) {
                        String codeStr = Integer.toBinaryString(codeCounter);
                        while (codeStr.length() < 32) {
                            codeStr = "0".concat(codeStr);
                        }
                        symbolTable.put(sKey, codeStr);
                        codeCounter++;
                    }
                }
            }
            System.out.println(codeCounter + " initial codes written to table");

            // LZW compression
            Boolean first = true;
            String s = "";
            String c ;
            int buffer = 0;
            int bitCounter = 0;

            // First write num of columns to file
            out.writeInt(columnCount);
            System.out.println(data[0].length + " columns");
            for (int j = 0; j < 20; j++) {
                for (int i = 0; i < data[j].length; i++) {
                    if (first) {
                        s = String.valueOf(data[j][0]);
                        c = String.valueOf(data[j][1]);
                        i = i + 1;
                        first = false;
                    } else {
                        c = String.valueOf(data[j][i]);
                    }
                    String t = s.concat(",");
                    t = t.concat(c);
                    if ((symbolTable.get(t)) != null) {
                        s = s.concat(",");
                        s = s.concat(c);
                    } else {
                        String symbol = symbolTable.get(s);
                        char[] ch = symbol.toCharArray();
                        for (int k = 0; k < ch.length; k++) {
                            if (bitCounter == 32) {
                                out.writeInt(buffer);
                                buffer = 0;
                                bitCounter = 0;
                                buffersWritten++;
                            }
                            buffer <<= 1;
                            if (ch[k] == '0') {
                                buffer |= 0;
                            } else {
                                buffer |= 1;
                            }
                            bitCounter++;
                        }
                        String codeStr = Integer.toBinaryString(codeCounter);
                        while (codeStr.length() < 32) {
                            codeStr = "0".concat(codeStr);
                        }
                        String n = s.concat(",");
                        n = n.concat(c);
                        symbolTable.put(s.concat(",").concat(c), codeStr);
                        s = c;
                        codeCounter++;

                    }
                }

            }

            // finally, output s
            String symbol = symbolTable.get(s);
            System.out.println("Finally outputting s: " + symbol);
            char[] ch = symbol.toCharArray();
            for (int k = 0; k < ch.length; k++) {
                if (bitCounter == 32) {
                    out.writeInt(buffer);
                    buffer = 0;
                    bitCounter = 0;
                    buffersWritten++;
                }
                buffer <<= 1;
                if (ch[k] == '0') {
                    buffer |= 0;
                } else {
                    buffer |= 1;
                }
                bitCounter++;
            }

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

        System.out.println("Buffers written: " + buffersWritten);

        return symbolTable;
    }

    public static double DO4(HashMap<String, String> symbolTable, int r) 
    {
        int buffer;
        int bitMask = 1 << 31;
        int columnCount = 0;
        int bitCounter = 0;
        String finalResult = "";

        try {
            DataInputStream dis = new DataInputStream(new FileInputStream("encode"));
            columnCount = dis.readInt();
            while (true) {
                String symbol = "";
                buffer = dis.readInt();
                for (int i = 0; i < 32; i++) {
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

                //	System.out.println("Looking up " + symbol);				
                for (String key : symbolTable.keySet()) {

                    if ((symbolTable.get(key)).equals(symbol)) {
                        //	System.out.println("Found " + symbol + " at " + key);
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
            System.out.println("Value list is: " + valueList.length);

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
    
    //Helper methods for numerical approximation for normal cdf with mean 0 and std dev = 0.33
    
    public static double cdf(double x)
    {
        double temp = x / 0.33;
        return Math.pow(1+Math.exp(0.0054-1.6101*temp-0.0674*Math.pow(temp, 3)),-1);
    }

    //Helper method for scaling and interpolating RGB values
    public static int getRGBColor(int c1, int c2, double in)
    {
        int c1r, c1b, c1g;
        int c2r, c2b, c2g;
        int c12r, c12b, c12g;
        int result;
        
        //unpack components in c1/c2
        c1b = c1&0x0000ff;
        c2b = c2&0x0000ff;
        c1g = (c1>>8)&0x0000ff;
        c2g = (c2 >> 8)&0x0000ff;
        c1r = (c1>>16)&0x0000ff;
        c2r = (c2 >> 16)&0x0000ff;
        
        //calculate direction vector
        c12r = c2r-c1r;
        c12g = c2g-c1g;
        c12b = c2b-c1b;
        
        //scale, repack, and return
        result = (int) (c1r+in*c12r);
        result = (result << 8) | (int) (c1g+in*c12g);
        result = (result << 8) | (int) (c1b+in*c12b); 
        
        return result;
    }

    //Tree Helper Functions
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
                parent.symbol = Integer.toBinaryString(parentSymbol);
                a.parent = parent;
                b.parent = parent;
                parent.left = b;
                parent.right = a;
                parentSymbol++;
                q.add(parent);
            }          
            return q.poll();
        }
    }

    public static Comparator<Node> nodeComparator = new Comparator<Node>() {
        @Override
        public int compare(Node n1, Node n2) {
            return (int) (n1.value - n2.value);
        }
    };

    // </editor-fold>
}
