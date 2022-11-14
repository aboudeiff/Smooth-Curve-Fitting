package smoothCurve;



import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Scanner;

public class Main {
	
	public static void main(String [] args) throws FileNotFoundException  {
	
		Scanner input =new Scanner(new FileReader("D:\\Software Engineering\\9. Ninth term\\Soft and Genetic\\Assignment2\\curve_fitting_input.txt"));
		int testCase=input.nextInt();
		
		// to perform number of test cases that entered
		for(int i=1;i<=testCase;i++)
		{
			int pointsNum=input.nextInt();
			int degree=input.nextInt();
			
			//array for points
			Point[] points=new Point[pointsNum];
			
			//insert points in array
			for(int j=0;j<pointsNum;j++)
			{
				Point point=new Point();
				point.x=input.nextDouble();
				point.y=input.nextDouble();
				points[j]=point;
			}
			
			System.out.println("The Test Case Number: "+i);
			
			//call constructor
			Chromosome chromosome=new Chromosome(degree);
			chromosome=Perform.run(degree, points);
		}
		input.close();
	}
}
