package smoothCurve;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class Perform {

	private static int populationSize = 4;
	private static int maxGeneration = 1;
	private static int selectionSize = 2;
	private static double lBound = -10.0;
	private static double uBound = 10.0;

	public static Chromosome run(int degree, Point[] points) {

		// Generate Random Population
		Chromosome[] population = makePopulation(populationSize, degree);

		for (int t = 0; t < maxGeneration; t++) {
			// Evaluate Fitness Function
			calculateFitness(population, points);
			// Perform Selection
			Chromosome[] selectedChromosomes = selection(selectionSize, 3, population);

			// Perform CrossOver
			Chromosome[] offsprings = crossOver(selectedChromosomes, points);

			// Perform Mutation
			mutation(offsprings, points, t, degree);

			// Perform Replacement
			replacement(population, points, offsprings);
		}

		// get best chromosome
		double fitness = population[0].fitness;
		int index = 0;

		for (int k = 1; k < population.length; k++) {
			if (fitness < population[k].fitness) {
				index = k;
				fitness = population[k].fitness;
			}
		}
		return population[index];
	}

//-----------------------------------------------------------------------------------------------------------------------
	// generate population
	private static Chromosome[] makePopulation(int size, int degree) {

		Chromosome[] population = new Chromosome[size];

		for (int i = 0; i < size; ++i) {
			population[i] = new Chromosome(degree);
			population[i].makeChromosome(lBound, uBound);
		}

		return population;
	}

//------------------------------------------------------------------------------------------------------------------------
	// perform fitness to all chromosomes
	private static void calculateFitness(Chromosome[] population, Point[] points) {
		for (int i = 0; i < population.length; ++i) {
			population[i].SingleChromosomeFitness(points);
		}
	}

//------------------------------------------------------------------------------------------------------------------------
	// perform selection ,
	// selection size = wanted parents , chromosomeTournamentSize = 3 (selected)
	// chromo for tournamnet)
	private static Chromosome[] selection(int selectionSize, int chromosomeTournamentSize, Chromosome[] population) {

		Chromosome[] selectedChromosomes = new Chromosome[selectionSize];// to fill the it with parents

		for (int i = 0; i < selectionSize; i++)// tournament selection loop ala hasab el parents
		{
			Random random = new Random();
			Chromosome[] chromosomeTournament = new Chromosome[chromosomeTournamentSize];// array for chromosome enter
																							// tour

			for (int j = 0; j < chromosomeTournamentSize; j++)// ala hasab K
				// take chromosomes from pop and enter tour
				chromosomeTournament[j] = population[random.nextInt(population.length)];

			int bestChromosome = 0;
			for (int j = 0; j < chromosomeTournamentSize; j++)// loop on chromosome that are in tournament ala hasab
																// fitness
			{
				if (chromosomeTournament[j].fitness > chromosomeTournament[bestChromosome].fitness) {
					bestChromosome = j;
				}
			}

			selectedChromosomes[i] = chromosomeTournament[bestChromosome];
		}

		System.out.println("----------------------After performing selection-------------------------");
		for (int i = 0; i < selectedChromosomes.length; i++) {
			System.out.println("Selected Chromosome: " + (i + 1));
			selectedChromosomes[i].print();
		}
		return selectedChromosomes;
	}

// ------------------------------------------------------------------------------------------------------------------------
	// selected chromosome as a parameter
	private static Chromosome[] crossOver(Chromosome[] selectedChromosomes, Point[] points) {

		// copied in offspring
		Chromosome[] offsprings = Arrays.copyOf(selectedChromosomes, selectedChromosomes.length);
		int point1 = 0,point2 = 0;

		// loop to perform crossover on all chromosome
		for (int i = 0; i < selectedChromosomes.length; i += 2) {

			Random r1 = new Random();
			// two points random
			point1 = r1.nextInt(selectedChromosomes[i].chromosome.length);
			point2 = r1.nextInt(selectedChromosomes[i].chromosome.length);
			if (point2 < point1) {
				int tmp = point2;
				point2 = point1;
				point1 = tmp;
			}

			// to check cross over
			double r2 = r1.nextDouble();
			if (r2 < 0.6) {

				Chromosome tmp = new Chromosome(offsprings[i]);

				// offspring call the crossover function and change
				offsprings[i].twoPointcrossOver(point1, point2, offsprings[i + 1]);
				offsprings[i].SingleChromosomeFitness(points);

				offsprings[i + 1].twoPointcrossOver(point1, point2, tmp);
				offsprings[i + 1].SingleChromosomeFitness(points);

			}
		}
		System.out.println("----------------------After performing crossover-------------------------");
		System.out.println("Crossover Point 1: " + point1);
		System.out.println("Crossover Point 2: " + point2);

		for (int i = 0; i < selectedChromosomes.length; i++) {
			System.out.println("Offspring: " + (i + 1));
			offsprings[i].print();
		}
		return offsprings;
	}

//------------------------------------------------------------------------------------------------------------------------
	// offsprings from crossover
	private static void mutation(Chromosome[] offsprings, Point[] points, int currentGeneration, int degree) {

		Random random = new Random();
		boolean flag = true;
		double r1 = 0.0, checkMutation_rm = 0.0, b = 0.5;
		double y, lowerBound_delta, upperBound_delta, mutation;

		// to fill new offsprings
		for (int i = 0; i < offsprings.length; i++) {

			// to loop on each gene in offspring and check mutation
			for (int j = 0; j < (degree + 1); j++) {

				checkMutation_rm = random.nextDouble();

				// check mutation with pm
				if (checkMutation_rm <= 0.1) {
					r1 = random.nextDouble();
					lowerBound_delta = offsprings[i].chromosome[j] - lBound;
					upperBound_delta = uBound - offsprings[i].chromosome[j];

					// check r1 with respect to r1
					if (r1 > 0.5) {
						y = upperBound_delta;
						flag = true;
					} else {
						y = lowerBound_delta;
						flag = false;
					}

					double r = random.nextDouble();
					// delta(t,y)
					mutation = y * (1 - Math.pow(r, Math.pow(1 - (currentGeneration / maxGeneration), b)));

					if (flag) {
						offsprings[i].chromosome[j] += mutation;
					} else {
						offsprings[i].chromosome[j] -= mutation;
					}

				}
			}

		}

		System.out.println("----------------------After performing Mutation-------------------------");
		for (int i = 0; i < offsprings.length; i++) {
			offsprings[i].SingleChromosomeFitness(points);
			System.out.println("Offspring: " + (i + 1));
			offsprings[i].print();
		}

	}

//------------------------------------------------------------------------------------------------------------------------
	private static ArrayList<Integer> findWorstChromosomes(Chromosome[] population) {
		ArrayList<Chromosome> staticTmp = new ArrayList<Chromosome>();
		for (Chromosome ch : population)
			staticTmp.add(ch);

		ArrayList<Chromosome> tmp = new ArrayList<Chromosome>();
		for (Chromosome ch : population)
			tmp.add(ch);

		ArrayList<Integer> worst = new ArrayList<Integer>();
		for (int i = 0; i < selectionSize; i++) {
			Chromosome worstChromosome = tmp.get(0);
			for (int j = 0; j < tmp.size(); j++) {
				if (tmp.get(j).fitness < worstChromosome.fitness)
					worstChromosome = tmp.get(j);
			}
			worst.add(staticTmp.indexOf(worstChromosome));
			tmp.remove(worstChromosome);
		}
		return worst;
	}

	private static void replacement(Chromosome[] population, Point[] points, Chromosome[] offsprings) {

		ArrayList<Integer> worstChromosomesinPopulation = findWorstChromosomes(population);
		ArrayList<Integer> worstChromosomesinOffsprings = findWorstChromosomes(offsprings);

		for (int i = 0; i < worstChromosomesinPopulation.size(); i++) {
			int selected = 0;
			for (int j = 0; j < worstChromosomesinOffsprings.size(); j++) {

				if (population[worstChromosomesinPopulation.get(i)].fitness < offsprings[worstChromosomesinOffsprings.get(j)].fitness) {
					population[worstChromosomesinPopulation.get(i)] = offsprings[worstChromosomesinOffsprings.get(j)];
					selected = j;
				}
			}
			worstChromosomesinOffsprings.remove(selected);

		}

		System.out.println("----------------------After performing Replacement-------------------------");
		System.out.println("Replaced Population ");

		for (int i = 0; i < population.length; i++) {
			population[i].SingleChromosomeFitness(points);
			population[i].print();
		}

	}

}
