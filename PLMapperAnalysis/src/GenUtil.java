import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class GenUtil {
    static public final DecimalFormat zeroDigit = new DecimalFormat("#0");
    static public final DecimalFormat oneDigit = new DecimalFormat("#0.0");
    static public final DecimalFormat twoDigit = new DecimalFormat("#0.00");
    static public final double DEGREES=57.29577951;   //   180.0/PI (radian -> degree)
    static public final double RAD = 0.01745329252; // PI/180.0 (degree -> radian)
    
	static public int coordinateToIndex(int coordinate, int maxRadius)
	{
		int index = coordinate + maxRadius;
		if((index < 0) ||(index > (2*maxRadius))) return -1;
		return index;
	}

	static public int indexToCoordinate(int index, int maxRadius)
	{
		return (index - maxRadius);
	}
    
	static public boolean isOutsideByIndex(int Xindex,int Yindex,int maxRadius) {
    	double cradius = Math.sqrt(1.0*GenUtil.indexToCoordinate(Xindex,maxRadius)*GenUtil.indexToCoordinate(Xindex,maxRadius) 
    			+ 1.0*GenUtil.indexToCoordinate(Yindex,maxRadius)*GenUtil.indexToCoordinate(Yindex,maxRadius));
    	
    	if(cradius > maxRadius) return true;
    	else return false;
	}
	
	static public boolean isOutsideByCoordinate(int Xcoord,int Ycoord,int maxRadius) {
    	double cradius = Math.sqrt(1.0* Xcoord* Xcoord + 1.0*Ycoord * Ycoord);
    	
    	if(cradius > maxRadius) return true;
    	else return false;
	}
	
	
    static public String [] strbufferToStringarr(StringBuffer buffer2,String termstr){

        StringBuffer buffer = reshapeStrbuffer(buffer2,termstr);
        int nelements = strbufferNElements(buffer);
        String [] ival = new String [nelements];

        int pos1 = 0,pos2 = 0, blen;
        int idx = 0;
        do{
            char posChar = buffer.charAt(pos1);
            while((posChar == ' ')||(posChar == '\t')||(posChar == ',')||(posChar == '<')||(posChar == '>')
                    ||(posChar == '=')||(posChar == '\"')||(posChar == ':')){
                if(pos1 == (buffer.length()-1)){
                    String [] ival2 = new String[idx];
                    for(int i= 0;i<idx;i++)
                        ival2[i] = ival[i];
                    return ival2;
                }
                buffer.deleteCharAt(pos1);
                posChar = buffer.charAt(pos1);
            }
            pos2 = pos1 + 1;
            blen = buffer.length();

            while(pos2 < blen){
                posChar = buffer.charAt(pos2);
                if((posChar != ' ')&&(posChar != '\t')&&(posChar != ',')&&(posChar != '<')&&(posChar != '>')
                        &&(posChar != '=')&&(posChar != '\"')&&(posChar != ':')){
                    pos2++;
                }
                else break;
            }

            try{
                ival[idx] = buffer.substring(pos1,pos2);
                idx++;
            }catch(NumberFormatException e){}

            if(pos2 >= blen){
                    String [] ival2 = new String[idx];
                    for(int i= 0;i<idx;i++)
                        ival2[i] = ival[i];
                    return ival2;
            }
            pos1 = pos2;

        }while(idx < ival.length);

        return ival;
    }
    
	@SuppressWarnings({ "rawtypes", "unchecked" })
	static public double[][] readDataFile(String filename,int[] selected,int[] colidx,boolean addxcol){
//	    int[] selected = new int[] {RowsToSkipInitially,RowsToRead,RowsToSkipInBetween};
//	    int[] colidx = new int[] {xidx,yidx1,yidx2,...};
	        int ncols = colidx[colidx.length-1] + 1;
	        ArrayList[] al = new ArrayList[ncols];
	        for(int i = 0; i < ncols; i++){
	            al[i] = new ArrayList();
	        }

	        StringBuffer buffer = new StringBuffer("");
			File readfile = new File(filename);
	        if(!readfile.exists()){
	            saveToFile("20_Error.log",(filename + " does not exist!!!"), false);
	            System.exit(0);
	        }

	        try{
	            BufferedReader fbr = new BufferedReader(new FileReader(readfile));
	            boolean cont = true;
	            for(int i = 0; i < selected[0]; i++) fbr.readLine();
	            do{
	                for(int i = 0; i < selected[1]; i++){
	                    buffer.append(fbr.readLine());
	                    if(buffer.length() == 0){
	                        cont = false;
	                        break;
	                    }
	                    double[] dval = new double [ncols];
	                    double [] dval2 = strbufferToDoublearr(buffer,"$$");
	                    if(dval2.length == 0) { cont = false; break;}
	                    else cont = true;
	                    buffer.delete(0, buffer.length());
	                    for(int j = 0; j < ncols; j++){
	                        if(j < dval2.length) dval[j] = dval2[j];
	                        al[j].add(dval[j]);
	                    }
	                }
	                if (!cont) break;

	                for(int i = 0; i < selected[2]; i++) fbr.readLine();
	                buffer.delete(0, buffer.length());
	            }while(cont);
	            fbr.close();
	        }catch(IOException ee){ System.exit(0);}

	        int nrows = al[0].size();
	        double[][] xy = new double[ncols][nrows];
	        for(int i = 0; i < ncols; i++){
	            for(int j = 0; j < nrows; j++){
	                xy[i][j] = Double.valueOf(al[i].get(j).toString()).doubleValue();
	            }
	        }

	        int totncols;
	        if(addxcol){
	            totncols = colidx.length + 1;
	        }
	        else{
	            totncols = colidx.length;
	        }

	        double[][] drawxy = new double[totncols][nrows];

	        for(int i = 0; i < totncols; i++){
	            for(int j = 0; j < nrows; j++){
	                if(addxcol){
	                    if(i == 0)
//	                        drawxy[i][j] = (double)(j+1);
	                        drawxy[i][j] = (double)j;
	                    else
	                        drawxy[i][j] = xy[colidx[i-1]][j];
	                }
	                else{
	                    drawxy[i][j] = xy[colidx[i]][j];
	                }
	            }
	        }

	        return drawxy;
	    }
    
    static public StringBuffer reshapeStrbuffer(StringBuffer bufferOriginal, String termchar){
        StringBuffer buffer = new StringBuffer(bufferOriginal.toString());
        StringBuffer retbuffer = new StringBuffer("");

        int pos = buffer.indexOf(termchar);
        if(pos >= 0){
            String tempstr = buffer.substring(0,pos);
            retbuffer.append(tempstr);
            return retbuffer;
        }
        else
            return buffer;
    }

    static public int strbufferNElements(StringBuffer buffer0){
        StringBuffer buffer = new StringBuffer("");
        buffer.append(buffer0);
        int pos1 = 0,pos2 = 0, blen;
        int idx = 0;
        do{
            while((buffer.charAt(pos1) == ' ')||(buffer.charAt(pos1) == '\t')||(buffer.charAt(pos1) == ',')
                    ||(buffer.charAt(pos1) == '=')||(buffer.charAt(pos1) == '\"')||(buffer.charAt(pos1) == ':')){
                if(pos1 == (buffer.length()-1))
                    return idx;
                buffer.deleteCharAt(pos1);
            }
            pos2 = pos1 + 1;
            blen = buffer.length();
            while((pos2 < blen)&&(buffer.charAt(pos2) != ' ')&&(buffer.charAt(pos2) != '\t')&&(buffer.charAt(pos2) != ',')
                    &&(buffer.charAt(pos2) != '=')&&(buffer.charAt(pos2) != '\"')&&(buffer.charAt(pos2) != ':')){
                pos2++;
            }

            idx++;
            if(pos2 >= blen)
                return idx;

            pos1 = pos2;
        }while(true);
    }
    
    static public double [] strbufferToDoublearr(StringBuffer buffer2,String termstr){
        StringBuffer buffer = reshapeStrbuffer(buffer2,termstr);
        int nelements = strbufferNElements(buffer);
        double [] dval = new double [nelements];

        int pos1 = 0,pos2 = 0, blen;
        int idx = 0;
        do{
            char posChar = buffer.charAt(pos1);
            while((posChar == ' ')||(posChar == '\t')||(posChar == ',')||(posChar == '<')||(posChar == '>')
                    ||(posChar == '=')||(posChar == '\"')||(posChar == ':')){
                if(pos1 == (buffer.length()-1)){
                    double [] dval2 = new double[idx];
                    for(int i= 0;i<idx;i++)
                        dval2[i] = dval[i];
                    return dval2;
                }
                buffer.deleteCharAt(pos1);
                posChar = buffer.charAt(pos1);
            }
            pos2 = pos1 + 1;
            blen = buffer.length();

            while(pos2 < blen){
                posChar = buffer.charAt(pos2);
                if((posChar != ' ')&&(posChar != '\t')&&(posChar != ',')&&(posChar != '<')&&(posChar != '>')
                        &&(posChar != '=')&&(posChar != '\"')&&(posChar != ':')){
                    pos2++;
                }
                else break;
            }

            try{
                dval[idx] = Double.parseDouble(buffer.substring(pos1,pos2));
                idx++;
            }catch(NumberFormatException e){}

            if(pos2 >= blen){
                double [] dval2 = new double[idx];
                for(int i= 0;i<idx;i++)
                    dval2[i] = dval[i];
                return dval2;
            }
            pos1 = pos2;

        }while(idx < dval.length);

        return dval;
    }

    static public double[][] readRowData(String filename,int nrowstoskip,int ndataskip,boolean readxdata){

        double [][] dval = null;

        StringBuffer buffer_x = new StringBuffer("");
        StringBuffer buffer_y = new StringBuffer("");

        try{

            BufferedReader fbr = new BufferedReader(new FileReader(new File(filename)));

            for(int i = 0; i < nrowstoskip; i++) fbr.readLine();

            double [] dvaltmp_x,dvaltmp_y;
            int datasize;
            if(readxdata){
                buffer_x.append(fbr.readLine());
                buffer_y.append(fbr.readLine());
                if((buffer_x.length() == 0) || (buffer_y.length() == 0)){
                    datasize = 0;
                    dvaltmp_x = null;
                    dvaltmp_y = null;
                }
                else{
                    dvaltmp_x = strbufferToDoublearr(buffer_x,"$$");
                    dvaltmp_y = strbufferToDoublearr(buffer_y,"$$");
                    datasize = dvaltmp_y.length;
                }
            }
            else{
                buffer_y.append(fbr.readLine());
                if(buffer_y.length() == 0){
                    datasize = 0;
                    dvaltmp_x = null;
                    dvaltmp_y = null;
                }
                else{
                    dvaltmp_x = null;
                    dvaltmp_y = strbufferToDoublearr(buffer_y,"$$");
                    datasize = dvaltmp_y.length;
                }
            }
            fbr.close();


            if(datasize == 0) return null;

            dval = new double[2][datasize - ndataskip];

            for(int i = 0; i < (datasize - ndataskip); i++){
                if(readxdata) dval[0][i] = dvaltmp_x[ndataskip+i];
                else dval[0][i] = i;
    //            else dval[0][i] = i + 1;
                dval[1][i] = dvaltmp_y[ndataskip+i];
            }
        }catch(IOException ee){ 					
        	ee.printStackTrace();
        }

        return dval;
    }
    
    static public DataPlotSingle plotData(String header,double [][] xydata,String filename){
        return ( new DataPlotSingle(xydata,header,null,null,null,null,filename));
    }

//    static public DataPlotSingle plotData(String header,double [][] xydata
//            ,int [] plottype, int [] circlesize,double [] xrange,double[] yrange,String filename){
//        return ( new DataPlotSingle(xydata,header,plottype,circlesize,xrange,yrange,filename));
//    }

    static public int val_rgb (double min,double max,double val){
        if (val < min) val = min;
        if (val > max) val = max;
        double h = 240.0*(max - val)/(max - min);
        //Hue [0,360], 0 for Red, 120 for Green, 240 for Blue
        double s = 1.0; //Saturation [0,1], 1 for pure color
        double v = 1.0; // Brightness [0,1], 0 for black
        double r,g,b;

        if(s==0) { //achromatic(grey)
            r = g = b = v;
        }

	h /= 60;			// sector 0 to 5
	double i = Math.floor(h);
	double f = h - i;			// factorial part of h
	double p = v * ( 1 - s );
	double q = v * ( 1 - s * f );
	double t = v * ( 1 - s * ( 1 - f ) );
	switch( (int)(i +0.1) ) {
		case 0:
			r = v;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = v;
			b = p;
			break;
		case 2:
			r = p;
			g = v;
			b = t;
			break;
		case 3:
			r = p;
			g = q;
			b = v;
			break;
		case 4:
			r = t;
			g = p;
			b = v;
			break;
		default:		// case 5:
			r = v;
			g = p;
			b = q;
			break;
	}

        r = r * 255;
        g = g * 255;
        b = b * 255;

        int rgb = 65536 * (int)r + 256 * (int)g + (int)b;

        return rgb;
    }

    static public void saveToFile(String filename,String info,boolean add){
        try{
            PrintWriter fpw = new PrintWriter(new BufferedWriter((new FileWriter(new File(filename),add)), 1024));
            fpw.print(info);
            fpw.close();
        }catch(IOException ee){
        	ee.printStackTrace();
        }
    }

    static public double[] stat(double[] var){
        if(var.length < 1) return null;
        int ndata = var.length,mini=0,maxi=0;
        double avg=0.0,std=0.0,min=var[0],max=var[0],sum=0.0,sumsq=0.0;

        for(int i = 0 ; i < ndata; i++){
            if(var[i] < min){ min = var[i]; mini = i; }
            if(var[i] > max){ max = var[i]; maxi = i; }
            sum = sum + var[i];
            sumsq = sumsq + var[i] * var[i];
        }

        avg = sum/ndata;
        std = Math.sqrt((sumsq - ndata * avg * avg)/(ndata - 1));

        return (new double[]{avg,std,min,max,(ndata+0.1),(mini+0.1),(maxi+0.1)});
    }
    
    static public double[] stat(double xn,double [] stat0){
    	
    	if(stat0 == null){
	        return (new double []{xn,0.0,xn,xn,1.1});
    	}
    	else if(stat0[4] < 1){
	        return (new double []{xn,0.0,xn,xn,1.1});
    	}
    	else{
	    	double anm1 = stat0[0],snm1 = stat0[1],min0 = stat0[2],max0 = stat0[3];
	    	int n = (int)(stat0[4]) + 1;
	    	double an = xn/n + (n-1)*anm1/n;
	    	double sn2 = (xn*xn + snm1*snm1*(n-2) + (n-1)*anm1*anm1 - n*an*an)/(n-1);
	    	double min = Math.min(xn,min0);
	    	double max = Math.max(xn,max0);
	        return (new double[]{an,Math.sqrt(sn2),min,max,(n + 0.1)});
    	}
    }
    
    static public double[] stat2D(double[][] var,int[][] id){
    	// id == null -> includes all data
    	// id[i][j] > -1 to be included
        double avg=0.0,std=0.0,min=0.0,max=0.0,sum=0.0,sumsq=0.0;
        boolean includethis;

        int ndata=0;
        for(int i = 0;i<var.length;i++){
            for(int j=0; j<var[i].length;j++){
            	includethis= false;
                if(id == null) includethis = true;
                else if(id[i][j]>-1) includethis = true;
                
                if(includethis){
                	if(ndata == 0){
                		min = var[i][j];
                		max = var[i][j];
                	}
                	
                    if(var[i][j] < min) min = var[i][j];
                    if(var[i][j] > max) max = var[i][j];
                    sum = sum + var[i][j];
                    sumsq = sumsq + var[i][j] * var[i][j];
                    ndata++;
                }
            }
        }
        
        avg = sum/ndata;
        std = Math.sqrt((sumsq - ndata * avg * avg)/(ndata - 1));
        return (new double[]{avg,std,min,max,(ndata+0.1)});
    }
	
    static public void saveMatrix(double[][] matrix,int maxRadius,String filename) {
    	
		try
		{
			FileWriter writer = new FileWriter(filename);
			writer.write("Y\\X(mm)");
			for(int i = -maxRadius;i<=maxRadius;i++) {
				writer.write("," + i);
			}
			writer.write("\n");
			
			for(int j = 0;j < (2*maxRadius+1);j++) {
				writer.write((j - maxRadius) + "");
				for(int i = 0;i < (2*maxRadius+1);i++) {
					writer.write("," + twoDigit.format(matrix[i][j]));
				}				
				writer.write("\n");
			}
			
			writer.close();
		}
		catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
   static public double [][] getRadialPlotData(XYZ[] xyz,double[][][] allxy,double [] angle){
	   double[][] allXYdata = new double[xyz.length*2][];
	   for(int i=0;i<xyz.length;i++) {
		   double [][] xydata =  getSingleRadialPlotData(xyz[i],allxy[i],angle);
		   allXYdata[2*i] = xydata[0];
		   allXYdata[2*i+1] = xydata[1];
	   }
	   return allXYdata;
   }

   static public double [][] getTangentialPlotData(XYZ[] xyz,double[][][] allxy,double radius){
	   double[][] allXYdata = new double[xyz.length*2][];
	   for(int i=0;i<xyz.length;i++) {
		   double [][] xydata =  getSingleTangentialPlotData(xyz[i],allxy[i],radius);
		   allXYdata[2*i] = xydata[0];
		   allXYdata[2*i+1] = xydata[1];
	   }
	   return allXYdata;
   }

   // returns the remainder after division by modulus
   static public double Mod(double inval, double modulus)
   {
      return inval - (Math.floor(inval/modulus)*modulus);
   }

   static private int[] plotRadialTangentialXiYi(double radius,double angle,int [] xy0,int [] xysize){
   	int[] XiYi = new int[2];
   	double x = radius * Math.sin(angle*RAD);
   	double y = radius * Math.cos(angle*RAD);
   	XiYi[0] = roundDoubleValue(x - xy0[0]);
   	XiYi[1] = roundDoubleValue(y - xy0[1]);
   	if(XiYi[0] < 0) XiYi[0] = 0;
   	if(XiYi[1] < 0) XiYi[1] = 0;
   	if(XiYi[0] >= xysize[0]) XiYi[0] = xysize[0] - 1;
   	if(XiYi[1] >= xysize[1]) XiYi[1] = xysize[1] - 1;
   	return XiYi;
   }
   
   static private int roundDoubleValue(double x) {
	   if(x > 0) return (int)(x + 0.5);
	   else return (int)(x - 0.5);
   }

   // angle[0]: (-180,0)
   // angle[1]: [0,180]
   
    static public double [][] getSingleRadialPlotData(XYZ xyz,double [][] xydata,double[] angle){
    	int [] xx = xyz.xcoord;
    	int [] yy = xyz.ycoord;
    	
    	int [] xy0 = new int[]{xx[0],yy[0]};
    	int [] xysize = new int[]{(int)(xx[xx.length -1] * 2 + 1.001),(int)(yy[yy.length -1] *2 + 1.001)};
    	int [] XiYi = new int[2];
    	
    	int nrad = yy.length / 2 + 1;
    	double [][] radial1 = new double[2][nrad];
    	double [][] radial2 = new double[2][nrad];
        
    	for(int i = 0; i < nrad; i++){
    		radial1[0][i] = i;
    		XiYi = plotRadialTangentialXiYi(radial1[0][i],angle[0],xy0,xysize);
    		radial1[1][i] = xydata[XiYi[0]][XiYi[1]];
    	}
    	
    	for(int i = 0; i < nrad; i++){
    		radial2[0][i] = i;
    		XiYi = plotRadialTangentialXiYi(radial2[0][i],angle[1],xy0,xysize);
    		radial2[1][i] = xydata[XiYi[0]][XiYi[1]];
    	}
        
        int i1 = radial1[0].length,i2 = radial2[0].length;
        
        double [][] radial = new double[2][i1+i2];
        
        for(int i = 0; i < i1;i++){
        	radial[0][i] = -radial1[0][i1 - i -1];
        	radial[1][i] = radial1[1][i1 - i -1];
        }
        
        for(int i = 0; i < i2; i++){
        	radial[0][i1+i] = radial2[0][i];
        	radial[1][i1+i] = radial2[1][i];
        }
        
        return radial;
    }
    
    static public double [][] getSingleTangentialPlotData(XYZ xyz,double [][] xydata,double radius){ //String radialTangentialProfilePrefix,
    	int [] xx = xyz.xcoord;
    	int [] yy = xyz.ycoord;
    	
    	int [] xy0 = new int[]{xx[0],yy[0]};
    	int [] xysize = new int[]{(int)(xx[xx.length -1] * 2 + 1.001),(int)(yy[yy.length -1] *2 + 1.001)};
    	int [] XiYi = new int[2];
    	
    	int ntan = 360;
    	double [][] tangential = new double[2][ntan];
    	for(int i = 0; i < ntan; i++){
    		tangential[0][i] = i;
    		XiYi = plotRadialTangentialXiYi(radius,tangential[0][i],xy0,xysize);
    		tangential[1][i] = xydata[XiYi[0]][XiYi[1]];
    	}
    	
        return tangential;
    }
    
}

class Filename{
	String fullpath = null;
	String directory = null;
	String withExt = null;
	String withoutExt = null;
	String type = null; // string between last underscore and extension. If filename is "Test_abc.csv" then type is "abc" 
	
	public Filename(String fullpath) {
		if((new File(fullpath)).exists()) {
			this.fullpath = fullpath;
			this.directory = (new File(fullpath)).getParent();
	    	this.withExt = fullpath.substring(directory.length()+1,fullpath.length());
	    	this.withoutExt = fullpath.substring(directory.length()+1,fullpath.lastIndexOf("."));
	    	this.type = withoutExt.substring(withoutExt.lastIndexOf("_")+1);
		}
	}

    public String withoutExtPlus(String addition) {
    	if (withoutExt == null) return null;
    	return (directory + "\\" + withoutExt + addition);
    }
}

class XYZ{
	public int maxRadius;
	public int matrixSize;
	public int nxy;
	public int [] x;
	public int [] y;
	public double [] z;
	public int [] xcoord;
	public int [] ycoord;
	
	public XYZ(int nxy,int maxRadius) {
		this.nxy = nxy;
		this.maxRadius = maxRadius;
		this.matrixSize = 2 * maxRadius + 1;
		x = new int [nxy];
		y = new int [nxy];
		z = new double [nxy];
		
		xcoord = new int [matrixSize];
		ycoord = new int [matrixSize];
		
		for(int i= 0;i<matrixSize;i++) {
			xcoord[i] = GenUtil.indexToCoordinate(i, maxRadius);
			ycoord[i] = GenUtil.indexToCoordinate(i, maxRadius);
		}
	}
	
	public double [][] calculateMatrix(){
        double [][] matrix = new double[matrixSize][matrixSize];
        int [][] nz = new int[matrixSize][matrixSize];
        int i,j;

        for (i =0 ;i < matrixSize;i++){
            for (j=0;j < matrixSize;j++){
                nz[i][j] = 0;
                matrix[i][j] = 0.0;
            }
        }

        for (int k = 0; k < x.length; k++) {
        	i = GenUtil.coordinateToIndex(x[k],maxRadius);
        	j = GenUtil.coordinateToIndex(y[k],maxRadius);
        	
        	nz[i][j] = nz[i][j] + 1;
            matrix[i][j] = (matrix[i][j]*(nz[i][j]-1.0) + z[k])/nz[i][j];
        }
        
        for (i =0 ;i < matrixSize;i++){
            for (j=0;j < matrixSize;j++){
            	if(nz[i][j] > 0) continue;
            	if(GenUtil.isOutsideByIndex(i,j,maxRadius)) continue;
            	
            	int n = 1,ii,jj,nsub;
            	do {
            		nsub = 0;
            		for(ii = i-n; ii <= i+n;ii++) {
            			if((ii < 0)||(ii >= matrixSize)) continue;
            			
            			jj = j+n;
            			if(jj < matrixSize) {
	            			if(nz[ii][jj]>0) {
	            				nsub++;
	            				matrix[i][j] = (matrix[i][j]*(nsub-1.0) + matrix[ii][jj])/nsub;
	            			}
            			}
            			
            			jj = j-n;
            			if(jj > 0) {
	            			if(nz[ii][jj]>0) {
	            				nsub++;
	            				matrix[i][j] = (matrix[i][j]*(nsub-1.0) + matrix[ii][jj])/nsub;
	            				
	            			}
            			}
            		}
            		
            		for(jj = j - n + 1; jj < j + n;jj++) {
            			if((jj < 0)||(jj >= matrixSize)) continue;
            			
            			ii = i+n;
            			if(ii < matrixSize) {
	            			if(nz[ii][jj]>0) {
	            				nsub++;
	            				matrix[i][j] = (matrix[i][j]*(nsub-1.0) + matrix[ii][jj])/nsub;
	            			}
            			}
            			
            			ii = i-n;
            			if(ii > 0) {
	            			if(nz[ii][jj]>0) {
	            				nsub++;
	            				matrix[i][j] = (matrix[i][j]*(nsub-1.0) + matrix[ii][jj])/nsub;
	            				
	            			}
            			}
            		}
            		
            		n++;
            	}while(nsub == 0);
            }
        }
		return matrix;
	}
}
