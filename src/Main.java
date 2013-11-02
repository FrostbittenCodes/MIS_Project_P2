
//Quick and dirty java implementation of phase one for CSE 408 project 1


import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import javax.imageio.ImageIO;

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
        //Quantize(fx, 3, "output.csv", 8);
        
        //Visualizer
        //VisualizeNoise(fx, "output.csv");
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
    
    //Currently doesn't calculate noise, but outputs the appropriate visual stuff
    public static void VisualizeNoise(String f1, String f2)
    // <editor-fold>
    {
        BufferedReader input1;
        BufferedReader input2;
        String outputfile = "output.bmp";
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
    
    // </editor-fold>
}
