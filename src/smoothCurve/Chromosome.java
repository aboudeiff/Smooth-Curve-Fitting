package smoothCurve;

public class Chromosome {
	public double chromosome[]; // the array of chromosome
	public double fitness;// to use it in

	// constructor
	public Chromosome(int D) {
		chromosome = new double[D + 1]; // initialized the size of chorosome array

	}

	// copy constructor
	public Chromosome(Chromosome c) {
		fitness = c.fitness;
		chromosome = new double[c.chromosome.length];

		// fill the chromosome
		for (int i = 0; i < chromosome.length; ++i) {
			chromosome[i] = c.chromosome[i];
		}
	}
	
//------------------------------------------------------------------------------------------------------------------------
	
	double randomWithRange(double min, double max) {
		double range = (max - min) + 1;
		return (double) (Math.random() * range) + min;
	}
	
	
//------------------------------------------------------------------------------------------------------------------------
	
	// function to make chromosome
	public void makeChromosome(double lBound, double hBound) {
		for (int i = 0; i < chromosome.length; ++i) {
			Double rand = randomWithRange(lBound, hBound);
			chromosome[i] = rand;
		}
	}

// ------------------------------------------------------------------------------------------------------------------------
	/*
	 * public void generateRandomGenes(double low,double high) { for(int
	 * i=0;i<genes.length;++i) { Random random = new Random();
	 * genes[i]=(high-low)*random.nextDouble()+low;// why we use (high-low)? } }
	 */
	/*
	 * // [ 1/points sub form 1 to points (y_calc - y_actual)^2 ] public void
	 * perfromSingleChromosomeFitness(testt.Point[] points) {
	 * 
	 * fitness = 0.0; double totalError = 0.0; int counter = 1; // to loop on point
	 * (x,y) for (int i = 0; i < points.length; ++i) { double error = 0; double
	 * equation1 = 0; double equation2 = 0; double equation3 = 0; double equation4 =
	 * 0; double equation5 = 0;
	 * 
	 * // to loop on chromosome for (int j = 0; j <= chromosome.length; ++j) {
	 * //error =
	 * Math.pow(chromosome[j]+chromosome[j+1]*points[i].x+chromosome[j+2]*Math.pow(
	 * points[i].y, 2.0)-points[i].y, 2.0); equation1 = chromosome[j]*points[i].x;
	 * equation2 = chromosome[j]*Math.pow(points[i].y , 2.0); equation3 =
	 * chromosome[j]+ equation1 + equation2; equation4 = equation3 - points[i].y;
	 * equation5 = Math.pow(equation4, 2.0); error = equation5; totalError = error;
	 * } } totalError /= points.length; // 1/points fitness = 1.0 / totalError; //
	 * like the lap }
	 */
//------------------------------------------------------------------------------------------------------------------------
	
	
	
	// (1/N Σ (y_calc. – y_actual)^2)
	public void SingleChromosomeFitness(Point[] points) {

		fitness = 0.0;
		for (int i = 0; i < points.length; ++i) { // loop on points
			int x_axes = 1;// x^0
			double y_calc = 0;

			for (int j = 0; j < chromosome.length; ++j) {// loop on chromosome
				y_calc = y_calc + x_axes * chromosome[j];// Error at (1,5) = ((1.95 * 1^0+ 8.16 * 1^1 + -2 * 1^2) – 5)^2
															// = 9.67
				x_axes *= points[i].x;// to increase the power
			}
			// summation of error at each node added
			fitness = fitness + Math.pow(y_calc - points[i].y, 2.0); // (y_calc - y-actual)^2 //actual from text file
		}
		fitness = fitness / points.length; // total error / points
		fitness = 1.0 / fitness;
	}
	
	
//------------------------------------------------------------------------------------------------------------------------
	
	
	// 56789 no change now calling this function
	// and change the first chromosome
	public void twoPointcrossOver(int point1, int point2, Chromosome c) {

		for (int i = 0; i < c.chromosome.length; ++i) {// 01-23-4 56-78-9 //01784 56239
			if (i >= point1 && i < point2) {// if true switch
				// original ,taba3 call //as parameter
				this.chromosome[i] = c.chromosome[i];
			}
		}

	}
//------------------------------------------------------------------------------------------------------------------------
	
	
	public void print() {

		System.out.println();
		System.out.println("Chromosomes: ");
		for (int i = 0; i < chromosome.length; i++) {
			System.out.print("C" + (i + 1) + "		" + chromosome[i] + "\n");
		}
		System.out.println();
		System.out.println("Fitness: " + fitness + "\n");

	}

}
