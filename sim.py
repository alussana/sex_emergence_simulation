#!/usr/bin/env python3

# TODO: environment can change the gene_effects (additive fitness values)

from random import random
from numpy.random import normal

CARRYING_CAPACITY = 10000
START_N_GENOMES = 1000
N_GENES = 10000
PLOIDY = 2
MUTATION_RATE = 5 / N_GENES
CROSSING_OVER = True
CROSSINGOVER_LENGTH_MEAN = 50
CROSSINGOVER_LENGTH_SD = 5
START_FRAC_SEX_GENOMES = 0.1

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
        self.genes = [[random() < 0.5 for i in range(n_genes)] for i in range(p)]
        self.meiosis = False

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

    '''
    def __init__(self, c: int, n: int, p: int, s: int, h_mean: int, h_sd: float, x: float):
        self.carrying_capacity = c
        self.n_genes = n
        self.ploidy = p
        self.start_n_genomes = s
        self.gene_effects = [random() - 0.5 for i in range(n)]
        self.population = [Genome(n, p) for i in range(s)]
        self.crossingover_length_mean = h_mean
        self.crossingover_length_sd = h_sd
        self.start_sex_frac = x
        
    def propagate(self):
        genomes_indexes = [i for i in range(len(self.population))]
        fitness_list = [self.compute_fitness(self.population[i]) for i in range(len(self.population))]
        propagation_list = [x for _, x in sorted(zip(fitness_list, genomes_indexes))]
        new_population = []
        while len(new_population) < self.carrying_capacity and len(propagation_list) > 0:
            # 0. determine if the first genome from the queue propagation_list can sexually reproduce
            if self.population[propagation_list[0]].meiosis:
                # 1. remove and return the first genome from the queue propagation_list
                # 2. identify next in the propagation_list that has meiosis, at index j
                # 3. if no j, move to next iteration
                # 4. if j, remove and return it from propagation_list
                parent_a = propagation_list.pop(0)
                j = 0
                while self.population[propagation_list[j]].meiosis == False and j < len(propagation_list):
                    j = j + 1
                if j >= len(propagation_list):
                    parent_b = propagation_list.pop(j)
                    # 5. determine crossing-over points
                    crossingover_points = []
                    x = 0
                    x = int(x + normal(loc=self.crossingover_length_mean,
                                       scale=self.crossingover_length_sd))
                    while x < self.n_genes:
                        crossingover_points.append(x)
                        x = int(x + normal(loc=self.crossingover_length_mean,
                                           scale=self.crossingover_length_sd))
                    # 6. generate 4 new genomes and append them in new_population
                    # TODO
            else:
                new_population.append(self.population[propagation_list[i]])
                new_population.append(self.population[propagation_list[i]])
            i = i + 1

    def start(self, n, x):
        ## TODO
        ## 1. initialize meiosis attribute in the population according to frequency x
        ## 2. loop over n generations
        ##      1. propagate
        ##      2. mutate

if __name__ == '__main__':
    main()
