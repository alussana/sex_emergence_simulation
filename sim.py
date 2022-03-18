#!/usr/bin/env python3

from copy import deepcopy
from random import random
from numpy.random import normal
from random import shuffle
import matplotlib.pyplot as plt

class Genome:
    '''
    A genome made by n genes and having p ploidy, where p = 2 is a diploid
    genome. In the current version of the software, only p = 2 is allowed.

    The genome is initialized randomly, each gene having a random dosage of
    the alleles.

    By default, a genome is initialized without the ability to sexually
    reproduce
    '''
    def __init__(self, n: int, p: int):
        self.n_genes = n
        self.ploidy = p
        self.genes = [[random() < 0.5 for i in range(n)] for i in range(p)]
        #phase_a = [random() < 0.5 for i in range(n)]
        #phase_b = phase_a.copy()
        #shuffle(phase_a)
        #self.genes = [phase_a, phase_b]
        self.meiosis = False

    def mutate(self, mutation_rate):
        for i in range(len(self.genes)):
            for j in range(len(self.genes[i])):
                if random() < mutation_rate:
                    self.genes[i][j] = not self.genes[i][j]

    def generate_crossingover_points(self, h_mean: float, h_sd: float):
        crossingover_points = []
        x = 0
        x = int(x + normal(loc=h_mean, scale=h_sd))
        while x < self.n_genes:
            crossingover_points.append(x)
            x = int(x + normal(loc=h_mean,
                               scale=h_sd))
        if len(crossingover_points) == 0:
            print(
            f'''
            W: no crossing-over points generated. Check the following variables:
                * mean number of genes between consecutive crossing over joints: {h_mean}
                * Std dev of the number of genes between consecutive crossing over joints: {h_sd}
                * Number of genes in the chromosome: {self.n_genes}
            ''')
        return(crossingover_points)

    def perform_crossingover(self, crossingover_points):
        new_genes = deepcopy(self.genes)
        crossingover_points.insert(0, 0)
        crossingover_points.append(self.n_genes)
        crossing_over_segments = list(zip(crossingover_points[:len(crossingover_points)], crossingover_points[1:]))
        for s in range(len(crossing_over_segments)):
            if s % 2 == 0:
                start = crossing_over_segments[s][0]
                end = crossing_over_segments[s][1]
                new_genes[0][start:end], new_genes[1][start:end] = new_genes[1][start:end], new_genes[0][start:end]
        return(new_genes)

    def generate_gametes(self, h_mean: float, h_sd: float):
        if self.meiosis == False:
            print('W: self.meiosis is False, but self.generate_gametes() was called anyway')
        self.crossingover_length_mean = h_mean
        self.crossingover_length_sd = h_sd
        # duplicate genome and determine crossing-over points separately for each one of the two copies
        crossingover_points_a = self.generate_crossingover_points(h_mean = h_mean, h_sd = h_sd)
        crossingover_points_b = self.generate_crossingover_points(h_mean = h_mean, h_sd = h_sd)
        # perform crossing-over separately for each one of the two copies
        new_genome_a = self.perform_crossingover(crossingover_points_a)
        new_genome_b = self.perform_crossingover(crossingover_points_b)
        # generate gametes (meiosis II)
        self.gametes = [new_genome_a[0], new_genome_a[1], new_genome_b[0], new_genome_b[1]]

class Simulation:
    '''
    An environment in which Genomes evolve.

    Parameters:

        * c: the carrying capacity, i.e. the maximum number of genomes that
          can be concurrently sustained by the environment
        * n: number of genes in the evolving Genomes
        * p: ploidy of the evolving Genomes, where p = 2 is a diploid genome.
             In the current version of the software, only p = 2 is allowed
        * s: number of genomes that exist at the start of the simulation
        * h_mean: mean number of genes between two crossing-over points
          NOTE: the real number will be slightly lower due to the truncation
          at the end of the chromosome
        * h_sd: standard deviation of the number of genes between two
          crossing-over points
        * x: starting fraction of genomes that can sexually reproduce
        * g: number of generation the simulation will run
        * m: mutation rate, i.e. the probability of an allele to mutate in the
             next generation
    
    The following assumptions hold:

        * The genome is constituted by only one chromosome
        * Genes are ordered, meaning that a distance between two genes A and B
          can be calculated in terms of how many genes separate A and B
        * In the current version of the software, the genes' effects are
          independent
        * In the current version of the software, the effect of allele dosage
          is additive
        * Each gene has two alleles. One allele is supposed to have no effect
          on the fitness, while the other can be beneficial or detrimental
        * A genome is initialized randomly, each gene having a random dosage of
          the alleles
        * There is no additional cost associated with sexual reproduction
    '''
    def __init__(self, c: int, n: int, p: int, s: int, h_mean: int, h_sd: float, x: float, g: int, m: float):
        self.carrying_capacity = c
        self.n_genes = n
        self.ploidy = p
        self.start_n_genomes = s
        self.gene_effects = [random() - 0.5 for i in range(n)]
        self.population = [Genome(n, p) for i in range(s)]
        self.crossingover_length_mean = h_mean
        self.crossingover_length_sd = h_sd
        self.start_sex_frac = x
        self.n_generations = g
        self.mutation_rate = m

    def compute_fitness(self, genome):
        f = 0
        for a in range(len(genome.genes)):
            for g in range(len(genome.genes[a])):
                f += self.gene_effects[g] * genome.genes[a][g]
        return(f)    

    def propagate(self):
        # TODO test propagate()
        genomes_indexes = [i for i in range(len(self.population))]
        fitness_list = [self.compute_fitness(self.population[i]) for i in genomes_indexes]
        propagation_list = [x for _, x in sorted(zip(fitness_list, genomes_indexes), reverse=True)]
        new_population = []
        while (len(new_population) < self.carrying_capacity) and (len(propagation_list) > 0):
            # 1. remove and return the first index from the queue propagation_list
            parent_a_idx = propagation_list.pop(0)
            parent_a = self.population[parent_a_idx]
            if parent_a.meiosis:
                # 2. identify next in the propagation_list that has meiosis, at index j
                j = 0
                while j < len(propagation_list):
                    if self.population[propagation_list[j]].meiosis == True:
                        break
                    else:
                        j = j + 1
                # 3. if no j, move to next iteration
                if j < len(propagation_list):
                    # 4. if j, remove and return it from propagation_list
                    parent_b_idx = propagation_list.pop(j)
                    parent_b = self.population[parent_b_idx]
                    # 5. perform meiosis in parent_a and parent_b
                    parent_a.generate_gametes(h_mean = self.crossingover_length_mean,
                                              h_sd = self.crossingover_length_sd)
                    parent_b.generate_gametes(h_mean = self.crossingover_length_mean,
                                              h_sd = self.crossingover_length_sd)
                    # 6. generate 4 new genomes and append them in new_population
                    shuffled_gametes_idx = [i for i in range(len(parent_a.gametes))]
                    shuffle(shuffled_gametes_idx)
                    for g in range(len(shuffled_gametes_idx)):
                      new_genes = [parent_a.gametes[g], parent_b.gametes[shuffled_gametes_idx[g]]]
                      new_genome = Genome(n=self.n_genes, p=self.ploidy)
                      new_genome.meiosis = True
                      new_genome.genes = new_genes
                      new_population.append(new_genome)                      
            else:
                new_population.append(parent_a) # duplicate the individual
                new_population.append(parent_a)
                new_population.append(parent_a)
                new_population.append(parent_a)
        self.population = new_population

    def mutate_genomes(self):
        for i in range(len(self.population)):
          self.population[i].mutate(self.mutation_rate)

    def update_history(self):
        n_sex_allele = 0
        for i in range(len(self.population)):
            if self.population[i].meiosis:
                n_sex_allele += 1
        self.history['n_individuals'].append(len(self.population))
        self.history['n_sex_allele'].append(n_sex_allele)

    def plot_counts(self, output_path):
        n_sex_allele = self.history['n_sex_allele']
        n_individuals = self.history['n_individuals']
        plt.bar([i for i in range(len(n_individuals))], n_individuals, label='mitosis')
        plt.bar([i for i in range(len(n_individuals))], n_sex_allele, label='meiosis')
        plt.xlabel('Generation')
        plt.ylabel('Individuals')
        plt.legend()
        plt.gca().spines['right'].set_color('none')
        plt.gca().spines['top'].set_color('none')
        plt.tight_layout()
        plt.savefig(output_path, dpi=200)

    def start(self):
        # initialize meiosis attribute in the population according to frequency x
        for i in range(self.start_n_genomes):
            if random() < self.start_sex_frac:
                self.population[i].meiosis = True
        # initialise history
        self.history = {'n_individuals': [],
                        'n_sex_allele': []}
        self.update_history()
        # run simulation over n generations
        for i in range(self.n_generations):
            self.propagate()
            self.mutate_genomes()
            self.update_history()

def main():

    CARRYING_CAPACITY = 2000
    START_N_GENOMES = 100
    N_GENES = 100
    PLOIDY = 2
    MUTATION_RATE = 5 / N_GENES
    CROSSING_OVER = True
    CROSSINGOVER_LENGTH_MEAN = 10
    CROSSINGOVER_LENGTH_SD = 1
    START_FRAC_SEX_GENOMES = 0.1
    N_GENERATIONS = 20

    sim = Simulation(
        c=CARRYING_CAPACITY,
        n=N_GENES,
        p=PLOIDY,
        s=START_N_GENOMES,
        h_mean=CROSSINGOVER_LENGTH_MEAN,
        h_sd=CROSSINGOVER_LENGTH_SD,
        x=START_FRAC_SEX_GENOMES,
        g=N_GENERATIONS,
        m=MUTATION_RATE
    )
    
    sim.start()
    sim.plot_counts('counts.png')

if __name__ == '__main__':
    main()
